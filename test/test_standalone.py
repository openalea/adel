"""ShootArchitecture without given tables: Adel's sample plant by default, plants generated from a few scalars
(tables.wheat_tables), each plant its own in one population; grown by ScenarioDevelopment (shoot_scene)."""
import numpy as np
import pytest

from openalea.adel.development import shoot_scene
from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.tables import sample_tables, wheat_tables
from openalea.metafspm import Planting

VARIABLES = ("organ_length", "visible_length", "organ_area", "inclination", "shape_class", "shape_mature_length")


def grow(scene, n_steps):
    for _ in range(n_steps):
        scene.run()
    return scene.populations[0].data_structure, scene.populations[0].instance.architecture


def test_the_sample_plant_by_default():
    axes, phytomers = sample_tables()
    default, _ = grow(shoot_scene(), 60)
    given, _ = grow(shoot_scene(axes_table=axes, phytomer_table=phytomers), 60)
    for name in VARIABLES + ("x2", "z2", "area"):
        np.testing.assert_array_equal(np.asarray(default.get(name)), np.asarray(given.get(name)), err_msg=name)


def test_the_sample_plant_with_echap_shapes():
    ds, architecture = grow(shoot_scene(leaves=LeafShapes.echap(xy_model="Soissons_byleafclass")), 60)
    blades = architecture.organ_names() == "blade"
    assert set(np.asarray(ds.get("shape_class"))[blades].tolist()) <= {1, 2}
    assert np.sum(ds.get("green_area")) > 0


def test_generated_tables():
    axes, phytomers = wheat_tables(leaf_number=11, tiller_number=2, longest_blade_length=0.2,
                                   longest_blade_width=0.015)
    assert axes["axe"].tolist() == ["MS", "T1", "T2"]
    assert axes["nf"].tolist() == [11., 8., 7.]
    main = phytomers.loc[phytomers["axe"] == "MS"]
    assert main["Ll"].max() * 0.01 == pytest.approx(0.2)
    assert main["Lw"].max() * 0.01 == pytest.approx(0.015)
    assert main["Ll"].idxmax() in main.index[7:9]                       # about 2.4 leaves below the flag leaf
    tiller = phytomers.loc[phytomers["axe"] == "T1"]
    np.testing.assert_allclose(tiller["Ll"].to_numpy(), 0.9 * main["Ll"].to_numpy()[:8])
    assert wheat_tables(leaf_number=11, tiller_number=2, longest_blade_length=0.2,
                        longest_blade_width=0.015)[0] is axes                # cached


def test_plants_of_their_own_in_one_population():
    plants = [Planting(), Planting(x=0.15), Planting(x=0.3)]
    own = [dict(leaf_number=9, tiller_number=0), dict(leaf_number=12, tiller_number=2, longest_blade_length=0.25),
           dict(leaf_number=12, tiller_number=2, longest_blade_length=0.25)]
    scene = shoot_scene(plants, plant_parameters=own, phyllochron=100.)
    ds, architecture = grow(scene, 150)
    axes, phytomers, numbers, offsets = architecture.tables()
    assert len(set(numbers.tolist())) == 2 and len(axes) == 1 + 3          # the same scalars, the same tables
    plant_of_axis = ds.map_to(np.arange(3, dtype=float), "Plant", "Axis", "broadcast").astype(int)
    names = architecture.axis_names()
    assert [sorted(names[plant_of_axis == p].tolist()) for p in range(3)] == [["MS"], ["MS", "T1", "T2"],
                                                                              ["MS", "T1", "T2"]]
    final = np.asarray(ds.get("final_leaf_number"))
    assert final[(plant_of_axis == 0) & (names == "MS")].tolist() == [9]
    assert final[(plant_of_axis == 1) & (names == "MS")].tolist() == [12]
    organ_plant = ds.map_to(np.arange(3, dtype=float), "Plant", "Organ", "broadcast").astype(int)
    longest = [np.asarray(ds.get("shape_mature_length"))[organ_plant == p].max() for p in range(3)]
    assert longest[1] == pytest.approx(0.25) and longest[2] == pytest.approx(0.25) and longest[0] < 0.2
    assert np.asarray(ds.get("phytomer_number"))[(plant_of_axis == 0)].tolist() == [9]


def test_sown_below_the_coleoptile_and_t1_at_its_top():
    """A plant is sown at 0.9 times its coleoptile's final length (its main stem's rank-1 sheath): its base below
    ground, the coleoptile's tip above; T1, borne on the main stem's phytomer 1, starts at the coleoptile's top, as
    T2, and keeps its place in the MPG."""
    axes, phytomers = sample_tables()
    coleoptile = phytomers.loc[(phytomers["axe"] == "MS") & (phytomers["n"] == 1), "Gl"].iloc[0] * 0.01
    ds, architecture = grow(shoot_scene(), 60)
    assert ds.mtg.property("z")[ds.entity_ids("Plant")[0]] == pytest.approx(-0.9 * coleoptile)
    names, kind = architecture.axis_names(), architecture._kinds()[0]
    node_axis, z1 = np.asarray(ds.owner("Axis")), np.asarray(ds.get("z1"))
    first = {name: z1[np.flatnonzero(node_axis == axis)[0]] for axis, name in enumerate(names)}
    assert first["MS"] == pytest.approx(-0.9 * coleoptile)
    assert first["T1"] == pytest.approx(first["T2"]) and first["T1"] > 0.
    assert np.asarray(ds.get("insertion_rank"))[names == "T1"].tolist() == [1]
