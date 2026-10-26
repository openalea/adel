"""The users' entry point (openalea.adel.simulation): its summary agrees with the plant model run by shoot_scene, the
thermal time follows the given temperature, and the tables and drawing hold the plant."""
import numpy as np
import pandas as pd
import pytest

from openalea.adel.development import shoot_scene
from openalea.adel.simulation import MAIN_PARAMETERS, _value_at, parameters, simulate

DAYS = 40


@pytest.fixture(scope="module")
def default():
    return simulate(days=DAYS)


def test_the_summary_is_the_plant_model_run_alone(default):
    alone = shoot_scene(start=0., thermal_time_step=15.)
    ds = alone.populations[0].data_structure
    alone.run()
    for _ in range(DAYS):
        alone.run()
    last = default.summary.iloc[-1]
    assert list(default.summary["day"]) == list(range(DAYS + 1))
    assert last["thermal_time"] == pytest.approx(15. * DAYS)
    leaf = np.asarray(ds.labels()) == "LeafElement"
    assert last["green_leaf_area"] == pytest.approx(float(np.sum(np.asarray(ds.get("green_area"))[leaf])), rel=1e-12)
    assert last["haun_stage"] == pytest.approx(float(np.asarray(ds.get("haun_stage"))[0]), rel=1e-12)
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


def test_draw(default):
    import matplotlib
    matplotlib.use("Agg")
    assert default.draw().collections
    default.draw(color_by="area")


def test_values_at_a_day():
    assert _value_at(3., 1.5) == 3.
    assert _value_at([0., 10., 20.], 1.5) == 15.
    assert _value_at(pd.Series([0., 4.], index=[0., 2.]), 0.5) == 1.
    assert _value_at(lambda day: 2 * day, 1.5) == 3.
