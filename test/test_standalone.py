"""ShootArchitecture without given tables: Adel's sample plant by default, plants generated from a few scalars
(tables.wheat_tables), each plant its own in one population; grown by ScenarioDevelopment (shoot_scene)."""
import os

import numpy as np
import pandas as pd
import pytest

from openalea.adel.development import shoot_scene
from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.tables import DATA, echap_shape_indices, sample_tables, wheat_tables, with_tillers
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


def test_cnwheat_shapes_are_echap_s_mercia_profiles_and_soissons_midribs():
    preset, given = LeafShapes.cnwheat(), LeafShapes.echap(sr_model="Mercia", xy_model="Soissons_byleafclass")
    assert len(preset.fitted) == len(given.fitted)
    for mine, theirs in zip(preset.fitted, given.fitted):
        for a, b in zip(mine, theirs):
            np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(preset.relative_area([1, 2], [0., 0.2], [1., 0.7]),
                                  given.relative_area([1, 2], [0., 0.2], [1., 0.7]))


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
    scene = shoot_scene(plants, plant_scenarios=own, phyllochron=100.)
    ds, architecture = grow(scene, 150)
    axes, phytomers, numbers, offsets = architecture.tables()
    assert len(set(numbers.tolist())) == 2 and len(axes) == 1 + 3          # the same scalars, the same tables
    plant_of_axis = ds.owner("Plant", of="Axis")
    names = architecture.axis_names()
    assert [sorted(names[plant_of_axis == p].tolist()) for p in range(3)] == [["MS"], ["MS", "T1", "T2"],
                                                                              ["MS", "T1", "T2"]]
    final = np.asarray(ds.get("final_leaf_number"))
    assert final[(plant_of_axis == 0) & (names == "MS")].tolist() == [9]
    assert final[(plant_of_axis == 1) & (names == "MS")].tolist() == [12]
    organ_plant = ds.owner("Plant", of="Organ")
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
    node_axis, z1 = np.asarray(ds.owner("Axis")), np.asarray(ds.get("display_z1"))      # where they are drawn
    first = {name: z1[np.flatnonzero(node_axis == axis)[0]] for axis, name in enumerate(names)}
    assert first["MS"] == pytest.approx(-0.9 * coleoptile)
    assert first["T1"] == pytest.approx(first["T2"]) and first["T1"] > 0.
    assert np.asarray(ds.get("insertion_rank"))[names == "T1"].tolist() == [1]



def _grown_as_alone(own, shared, n_steps=80):
    """Plants of one population with their own scenarios *own* (and *shared*) grow as each alone."""
    ds, _ = grow(shoot_scene([Planting(x=0.15 * k) for k in range(len(own))], plant_scenarios=own, **shared), n_steps)
    organ_plant = ds.owner("Plant", of="Organ")
    node_plant = np.asarray(ds.owner("Plant"))
    triangles = ds.triangles()
    plant_areas = ds.map_to(triangles.per_owner(triangles.areas(), ds.n_nodes()), "node", "Plant", "sum")
    for plant, scenario in enumerate(own):
        alone, _ = grow(shoot_scene(**dict(shared, **scenario)), n_steps)
        for name in ("inclination", "organ_area", "shape_mature_length", "diameter", "organ_length"):
            np.testing.assert_allclose(np.sort(np.asarray(ds.get(name))[organ_plant == plant]),
                                       np.sort(np.asarray(alone.get(name))), rtol=1e-12, err_msg=f"{plant} {name}")
        np.testing.assert_allclose(np.sort(np.asarray(ds.get("area"))[node_plant == plant]),
                                   np.sort(np.asarray(alone.get("area"))), rtol=1e-12)
        assert plant_areas[plant] == pytest.approx(np.sum(alone.triangles().areas()), rel=1e-9)


def test_each_plant_its_own_phyllochrons_and_length_unit():
    """Each plant's phyllochrons (its blades' age classes in a dynamic database) and its tables' length unit (here a
    plant's tables in mm), in one population."""
    axes, phytomers = sample_tables()
    in_mm = echap_shape_indices(axes, phytomers)                       # the shapes of echap's database
    in_mm = in_mm.copy()
    in_mm[["Ll", "Lw", "Gl", "Gd", "El", "Ed"]] *= 10.
    axes_in_mm = axes.copy()
    axes_in_mm["dredT"] *= 10.
    _grown_as_alone([dict(), dict(phyllochron=100., tiller_phyllochron=110.),
                     dict(axes_table=axes_in_mm, phytomer_table=in_mm, length_unit=0.001)],
                    dict(leaves=LeafShapes.echap(xy_model="Soissons_byleafclass")))


def test_each_plant_its_own_inclination_rule():
    """Each plant's inclination rule (its tillers inclined or not, its ramp's start and rate), in one population."""
    _grown_as_alone([dict(), dict(inclination_rate=10., inclination_start_haun_stage=2.)],
                    dict(tiller_inclination=1))
    _grown_as_alone([dict(tiller_inclination=0), dict(tiller_inclination=1, inclination_rate=10.)], dict())


def _keyed(scene):
    """Per element of a grown *scene*, sorted by its key (plant, axis name, phytomer rank, organ, index in organ):
    (the keys, each one's topological parent's key, {variable: values})."""
    ds = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture
    g = ds.structure
    rank = np.asarray(ds.get("phytomer_rank"))[np.asarray(ds.owner("Phytomer"))]
    keys = list(zip(np.asarray(ds.owner("Plant")).tolist(),
                    architecture.axis_names()[np.asarray(ds.owner("Axis"))].tolist(), rank.tolist(),
                    architecture.organ_names()[np.asarray(ds.owner("Organ"))].tolist(),
                    np.asarray(ds.get("element_index")).tolist()))
    nodes = ds.entity_ids("node").tolist()
    position = {vid: i for i, vid in enumerate(nodes)}
    order = sorted(range(len(keys)), key=lambda i: keys[i])
    parents = [keys[position[g.parent(nodes[i])]] if g.parent(nodes[i]) in position else None for i in order]
    names = ("length", "area", "green_area", "senesced_area", "rolled_part", "width", "x1", "y1", "z1", "x2", "y2",
             "z2", "frame_x", "frame_y", "frame_z", "heading_x", "heading_y", "heading_z", "up_x", "up_y", "up_z")
    values = {name: np.asarray(ds.get(name))[order] for name in names}
    for name in ("organ_length", "visible_length", "inclination", "age", "organ_area"):
        values[name] = np.asarray(ds.get(name))[np.asarray(ds.owner("Organ"))][order]
    return [keys[i] for i in order], parents, values


def _tiller_of_tiller_tables():
    """Adel's sample plant with T1.1 (on T1's phytomer 1, its first leaf at 420 °Cd)."""
    axes = pd.read_csv(os.path.join(DATA, "axeTable.csv"))
    tiller = axes.loc[axes["axe"] == "T3"].assign(axe="T1.1", emf1=420., nf=5., nf_end=5.)
    axes = pd.concat([axes, tiller], ignore_index=True).sort_values("emf1", kind="stable").reset_index(drop=True)
    return axes, with_tillers(axes, pd.read_csv(os.path.join(DATA, "phytoT.csv")))


@pytest.mark.parametrize("own_axis_top", [0, 1])
def test_one_large_step_as_small_steps(own_axis_top):
    """A call creates every axis, phytomer and element its thermal time asks for (ScenarioShoot.grow): one large
    step (from the start, or after a first call) reaches the structure of small steps, its elements' values the
    same; plants with tillers, a tiller of a tiller, inclined."""
    axes, phytomers = _tiller_of_tiller_tables()
    own = [dict(leaf_number=12, tiller_number=3, longest_blade_length=0.25), dict(axes_table=axes,
                                                                                  phytomer_table=phytomers)]

    def scene(**parameters):
        return shoot_scene([Planting(), Planting(x=0.2)], plant_scenarios=own, delays=[0., 30.],
                           tiller_inclination=1, own_axis_top=own_axis_top, **parameters)

    small = scene(start=0., thermal_time_step=10.)
    for _ in range(81):                                  # the first call starts the thermal time, 80 add 10 °Cd
        small.run()
    after_first = scene(start=0., thermal_time_step=800.)
    for _ in range(2):
        after_first.run()
    at_once = scene(start=800.)
    at_once.run()
    expected_keys, expected_parents, expected = _keyed(small)
    assert any(key[1] == "T1.1" for key in expected_keys) and len({key[1] for key in expected_keys}) == 6
    for large in (after_first, at_once):
        keys, parents, values = _keyed(large)
        assert keys == expected_keys and parents == expected_parents
        for name, array in expected.items():
            np.testing.assert_allclose(values[name], array, rtol=1e-12, atol=1e-15, err_msg=name)


def test_blade_elements_from_their_base_to_their_end():
    """A blade element's drawn base (display_x1, ...) and end (display_x2, ...) lie along its blade, on its triangles' extent (its
    rolled part's cylinder and its flat part's mesh, drawn from its frame's origin frame_x…: its blade's collar for
    the first one; within 0.1 mm, its mesh's tip triangle left out below the minimal area), its base the previous
    element's end, its end away from its frame's origin. Its topological ends (x1 ... z2) keep that chord, each
    element's base at its parent's end (metafspm's ds.check_connectivity)."""
    ds, architecture = grow(shoot_scene(leaves=LeafShapes.echap(xy_model="Soissons_byleafclass"), n_sect=3,
                                        tiller_inclination=1), 70)
    triangles = architecture.triangles()
    leaf = architecture._kinds()[0] == 2
    base = np.stack([np.asarray(ds.get(f"display_{name}")) for name in ("x1", "y1", "z1")], axis=1)
    end = np.stack([np.asarray(ds.get(f"display_{name}")) for name in ("x2", "y2", "z2")], axis=1)
    topological = [np.stack([np.asarray(ds.get(name + k)) for name in "xyz"], axis=1) for k in "12"]
    np.testing.assert_allclose((topological[1] - topological[0])[leaf], (end - base)[leaf], atol=1e-12)
    assert ds.check_connectivity(tolerance=1e-12).empty
    frame = np.stack([np.asarray(ds.get(name)) for name in ("frame_x", "frame_y", "frame_z")], axis=1)
    drawn = np.unique(triangles.owner[leaf[triangles.owner]])
    assert drawn.size > 20
    for node in drawn.tolist():
        vertices = triangles.vertices[triangles.owner == node].reshape(-1, 3)
        low, high = vertices.min(axis=0) - 1e-4, vertices.max(axis=0) + 1e-4
        for point in (base[node], end[node]):
            assert ((point >= low) & (point <= high)).all(), node
    organ, index = np.asarray(ds.owner("Organ")), np.asarray(ds.get("element_index"))
    element = {(organ[node], index[node]): node for node in drawn.tolist()}
    rolled = np.asarray(ds.get("rolled_part"))
    pairs = [(element[(organ[node], index[node] - 1)], node) for node in drawn.tolist()
             if (organ[node], index[node] - 1) in element and rolled[node] == 0]
    assert pairs
    previous, following = np.array(pairs).T
    np.testing.assert_allclose(base[following], end[previous], atol=1e-9)
    assert (np.linalg.norm(end[drawn] - frame[drawn], axis=1) > 0).all()
    assert end[drawn, 2].max() > frame[drawn, 2].max() + 0.01             # the blades' tips above their collars


def test_a_stand_grown_to_its_stage_before_the_first_step():
    """A call before a scene's first step grows the plants to their stage (their start, minus their delays): the
    scene's first step adds no degree-days (none passed since), the next ones the step before's (15 °C, a day)."""
    scene = shoot_scene([Planting(), Planting(x=0.2)], delays=[0., 30.], time_step=86400.,
                        forcings={"air_temperature": 15.}, start=600.)
    model = scene.populations[0].instance
    model.grow()
    ds = scene.populations[0].data_structure
    grown = {name: np.asarray(ds.get(name)).copy() for name in ("organ_length", "z2")}
    np.testing.assert_allclose(np.asarray(ds.get("thermal_time")), [600., 570.])
    scene.run()
    np.testing.assert_allclose(np.asarray(ds.get("thermal_time")), [600., 570.])
    for name, values in grown.items():
        np.testing.assert_array_equal(np.asarray(ds.get(name)), values, err_msg=name)
    scene.run()
    np.testing.assert_allclose(np.asarray(ds.get("thermal_time")), [615., 585.])


def test_the_triangles_carry_their_green_share():
    """Each triangle's "green" attribute: the share of its area in its element's green part (the first green_length
    from the element's base, along a stem's axis or a blade's midrib): an element's triangles' green area, over their
    area, is its green_area over its area (within its mesh's resolution); a wholly green element's 1, a senesced
    one's 0."""
    ds, architecture = grow(shoot_scene(leaves=LeafShapes.echap(xy_model="Soissons_byleafclass"), n_sect=3), 90)
    triangles = architecture.triangles()
    green = triangles.attributes["green"]
    assert green.min() >= 0. and green.max() <= 1.
    n = ds.n_nodes()
    area = triangles.per_owner(triangles.areas(), n)
    green_area = triangles.per_owner(triangles.areas() * green, n)
    expected_area, expected_green = np.asarray(ds.get("area")), np.asarray(ds.get("green_area"))
    drawn = (area > 0) & (expected_area > 0)
    partly = drawn & (expected_green > 0) & (expected_green < 0.999 * expected_area)
    assert partly.sum() >= 5 and (drawn & (expected_green == 0)).any()
    np.testing.assert_allclose(green_area[drawn] / area[drawn], expected_green[drawn] / expected_area[drawn],
                               atol=0.05)
