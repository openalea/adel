"""The users' entry point (openalea.adel.simulation): its summary agrees with the plant model run by shoot_scene, the
thermal time follows the given temperature, and the tables and drawing hold the plant."""
import numpy as np
import pandas as pd
import pytest

from openalea.adel.development import ScenarioShoot, shoot_scene
from openalea.adel.simulation import MAIN_PARAMETERS, parameters, simulate
from openalea.metafspm.coupling.composite_wrapper import CompositeModel

DAYS = 40


@pytest.fixture(scope="module")
def default():
    return simulate(days=DAYS)


def test_the_summary_is_the_plant_model_run_alone(default):
    alone = shoot_scene(start=0., thermal_time_step=15.)
    data_structure = alone.populations[0].data_structure
    alone.run()
    for _ in range(DAYS):
        alone.run()
    last = default.summary.iloc[-1]
    assert list(default.summary["day"]) == list(range(DAYS + 1))
    assert last["thermal_time"] == pytest.approx(15. * DAYS)
    leaf = np.asarray(data_structure.labels()) == "LeafElement"
    assert last["green_leaf_area"] == pytest.approx(float(np.sum(np.asarray(data_structure.get("green_area"))[leaf])), rel=1e-12)
    assert last["haun_stage"] == pytest.approx(float(np.asarray(data_structure.get("haun_stage"))[0]), rel=1e-12)
    assert default.summary["tillers"].min() == 0.


def test_tables(default):
    organs, axes = default.organs(), default.axes()
    assert set(organs["organ"]) == {"internode", "sheath", "blade"}
    blades = organs[organs["organ"] == "blade"]
    assert blades["area"].sum() == pytest.approx(
        default.summary.iloc[-1]["green_leaf_area"] + default.summary.iloc[-1]["senesced_leaf_area"], rel=1e-9)
    assert axes["axis"].iloc[0] == "MS" and len(axes) == 1 + default.summary.iloc[-1]["tillers"]


def test_temperature_drives_the_development():
    warm = simulate(days=20, temperature=20.)
    series = simulate(days=20, temperature=pd.Series([10., 10.], index=[0., 20.]), base_temperature=5.)
    assert warm.summary.iloc[-1]["thermal_time"] == pytest.approx(400.)
    assert series.summary.iloc[-1]["thermal_time"] == pytest.approx(100.)
    assert series.summary.iloc[-1]["haun_stage"] < warm.summary.iloc[-1]["haun_stage"]


def test_parameters():
    table = parameters({"phyllochron": 100.})
    assert tuple(table["name"][:len(MAIN_PARAMETERS)]) == MAIN_PARAMETERS
    assert table.loc[table["name"] == "phyllochron", "value"].item() == 100.
    assert table["name"].is_unique and "thermal_time_step" not in set(table["name"])
    slower = simulate({"phyllochron": 120.}, days=20)
    assert slower.summary.iloc[-1]["haun_stage"] < simulate(days=20).summary.iloc[-1]["haun_stage"]


def test_a_stand():
    stand = simulate(days=20, plants=3)
    assert stand.axes()["plant"].nunique() == 3


class NestedShoot(CompositeModel):
    """ScenarioShoot nested in another plant model (as tutorial 6's LitShoot): its architecture a nested part's."""
    initiators = (ScenarioShoot,)
    parts = {"shoot": ScenarioShoot}


def test_a_model_nesting_the_shoot():
    nested, alone = simulate(days=20, model=NestedShoot), simulate(days=20)
    assert nested.architecture is nested.scene.populations[0].instance.shoot.architecture
    pd.testing.assert_frame_equal(nested.summary, alone.summary)


def test_draw(default):
    import matplotlib
    matplotlib.use("Agg")
    assert default.draw().collections
    default.draw(color_by="area")




def test_the_run_configuration_is_the_scenes():
    """compute= (metafspm's Compute) is the Scene's run configuration."""
    from openalea.metafspm.compute import Compute
    compute = Compute(threads=1, memory="keep")
    assert simulate(days=1, compute=compute).scene.compute is compute


def _every_value(scene) -> dict:
    """Every numeric variable of every population, by population and name."""
    found = {}
    for population in scene.populations:
        data_structure = population.data_structure
        for name in data_structure.available_vars():
            values = np.asarray(data_structure.get(name))
            if values.dtype != object:
                found[f"{population.name} {name}"] = values.copy()
    return found


def _two_populations_at(processes, build, steps) -> dict:
    """*build*(Compute)'s scene run *steps* steps: every value of its populations."""
    from openalea.metafspm.compute import Compute
    scene = build(Compute(threads=1, processes=processes))
    try:
        for _ in range(steps):
            scene.run()
        assert type(scene.runner).__name__ == ("ProcessRunner" if processes > 1 else "SequentialRunner")
        return _every_value(scene)
    finally:
        scene.finish()


def _assert_same_bits(one: dict, two: dict) -> None:
    assert one.keys() == two.keys() and one
    for key, values in one.items():
        np.testing.assert_array_equal(two[key], values, err_msg=key)


def test_two_populations_in_two_processes_give_the_bits_of_one():
    """Two plants in two populations (planting(population=)): Compute(processes=2) runs each in its worker process,
    with the values of one process."""
    from openalea.metafspm.scene.population import Planting, planting
    from openalea.metafspm.scene.scene import Scene

    def build(compute):
        table = planting([Planting(), Planting(x=0.2)], ScenarioShoot, {"start": 0.}, [{}, {"phyllochron": 120.}],
                         population=["early", "late"])
        return Scene(table, time_step=86400., forcings={"air_temperature": lambda time: 15.}, compute=compute)
    _assert_same_bits(_two_populations_at(1, build, 20), _two_populations_at(2, build, 20))
