"""ScenarioDevelopment (its default timings the harness' scenario) against legacy_run.given_lengths: every axis'
phytomer
number and every organ's lengths and age, at each step, on two plants (one delayed). And the blades' dimensions as
ShootArchitecture's inputs: a development model's values replace the phytomer table's."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import StandPlant, given_lengths, harness_tables, start
from openalea.adel.comparison.port_run import axis_names, organ_names, plant_of, port_scene

PLANTS = [StandPlant(), StandPlant(x=10., delay=37.)]
THERMAL_TIME_STEP = 10.


def expected(scene, lengths):
    """{variable: values} from given_lengths' records (*lengths*: one per plant), as cnwgrass writes them: an organ
    not appeared yet keeps its values (0)."""
    data_structure = scene.populations[0].data_structure
    names = axis_names(scene)
    axis_plant = plant_of(data_structure, "Axis")
    organ_axis = data_structure.owner("Axis", of="Organ")
    organ_plant, ranks, organs = plant_of(data_structure, "Organ"), np.asarray(data_structure.get("organ_rank")), organ_names(data_structure)
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
    data_structure = scene.populations[0].data_structure
    created = 0
    for step in range(80):
        scene.run()
        lengths = [given_lengths(start(axes) + step * THERMAL_TIME_STEP - plant.delay, axes, phytomers)
                   for plant in PLANTS]
        for name, values in expected(scene, lengths).items():
            np.testing.assert_allclose(np.asarray(data_structure.get(name), dtype=float), values, rtol=0., atol=1e-10,
                                       err_msg=f"{name} at step {step}")
        created = len(data_structure.entity_ids("Organ"))
    assert created > 30                                # the tillers appeared, phytomers were created


def test_the_blades_dimensions_from_development():
    scene = port_scene("samples", thermal_time_step=THERMAL_TIME_STEP)
    data_structure = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture
    for _ in range(25):
        scene.run()
    blade = organ_names(data_structure) == "blade"
    area_before = np.asarray(data_structure.get("organ_area"), dtype=float).copy()
    assert area_before[blade].max() > 0
    given = {}
    for name, shape in (("final_blade_length", "shape_mature_length"), ("final_blade_width", "shape_max_width")):
        values = np.full(len(blade), np.nan)                           # what a development model gives
        values[blade] = np.asarray(data_structure.get(shape), dtype=float)[blade] * 1.25
        data_structure.set(name, values)
        given[shape] = values
    architecture()
    area_after = np.asarray(data_structure.get("organ_area"), dtype=float)
    grown = blade & (area_before > 0)
    assert np.all(area_after[grown] != area_before[grown])
    for shape, values in given.items():                                # the blades' shapes: the given dimensions
        np.testing.assert_array_equal(np.asarray(data_structure.get(shape), dtype=float)[blade], values[blade])
    architecture()
    np.testing.assert_array_equal(np.asarray(data_structure.get("shape_max_width"), dtype=float)[blade],
                                  given["shape_max_width"][blade])      # kept
