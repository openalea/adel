"""ScenarioDevelopment (its default timings the harness' scenario) against legacy_run.given_lengths: every axis' phytomer
number and every organ's lengths and age, at each step, on two plants (one delayed). And the blades' dimensions as
ShootArchitecture's inputs: a development model's values replace the phytomer table's."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import given_lengths, harness_tables, start
from openalea.adel.comparison.port_run import axis_names, organ_names, plant_of, port_scene
from openalea.adel.development import Planting

PLANTS = [Planting(), Planting(x=10., delay=37.)]
THERMAL_TIME_STEP = 10.


def expected(scene, lengths):
    """{variable: values} from given_lengths' records (*lengths*: one per plant), as cnwgrass writes them: an organ
    not appeared yet keeps its values (0)."""
    ds = scene.populations[0].data_structure
    names = axis_names(scene)
    axis_plant = plant_of(ds, "Axis")
    organ_axis = ds.map_to(np.arange(names.size, dtype=float), "Axis", "Organ", "broadcast").astype(int)
    organ_plant, ranks, organs = plant_of(ds, "Organ"), np.asarray(ds.get("organ_rank")), organ_names(ds)
    out = dict(axis_number=np.array([float(sum(len(axis) > 0 for axis in mine.values())) for mine in lengths]),
               phytomer_number=np.array([len(lengths[plant][name]) for plant, name in zip(axis_plant, names)],
                                        dtype=float))
    for variable, given in (("organ_length", "length"), ("visible_length", "visible_length"),
                            ("organ_senesced_length", "senesced_length"), ("age", "age")):
        values = np.zeros(organs.size)
        for i in range(organs.size):
            organ = lengths[organ_plant[i]][names[organ_axis[i]]].get(int(ranks[i]), {}).get(organs[i])
            if organ is not None and given in organ:
                values[i] = organ[given] * (0.01 if given != "age" else 1.)
        out[variable] = values
    return out


def test_the_harness_scenario():
    axes, phytomers = harness_tables("samples")
    scene = port_scene("samples", plants=PLANTS, thermal_time_step=THERMAL_TIME_STEP)
    ds = scene.populations[0].data_structure
    created = 0
    for step in range(80):
        scene.run()
        lengths = [given_lengths(start(axes) + step * THERMAL_TIME_STEP - plant.delay, axes, phytomers)
                   for plant in PLANTS]
        for name, values in expected(scene, lengths).items():
            np.testing.assert_allclose(np.asarray(ds.get(name), dtype=float), values, rtol=0., atol=1e-10,
                                       err_msg=f"{name} at step {step}")
        created = len(ds.entity_ids("Organ"))
    assert created > 30                                # the tillers appeared, phytomers were created


def test_the_blades_dimensions_from_development():
    scene = port_scene("samples", thermal_time_step=THERMAL_TIME_STEP)
    ds = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture
    for _ in range(25):
        scene.run()
    blade = organ_names(ds) == "blade"
    area_before = np.asarray(ds.get("organ_area"), dtype=float).copy()
    assert area_before[blade].max() > 0
    for name in ("shape_mature_length", "shape_max_width"):
        values = np.asarray(ds.get(name), dtype=float).copy()
        values[blade] *= 1.25
        ds.set(name, values)
    architecture()
    written = np.asarray(ds.get("shape_max_width"), dtype=float).copy()
    area_after = np.asarray(ds.get("organ_area"), dtype=float)
    grown = blade & (area_before > 0)
    assert np.all(area_after[grown] != area_before[grown])
    architecture()
    np.testing.assert_array_equal(np.asarray(ds.get("shape_max_width"), dtype=float), written)   # kept
