"""The tillers' inclinations by static Adel's rules (tiller_inclination=1) against a line by line transcription of
Adel.R's getIncT and axe_inclination (R does not run here): every stem organ's inclination, every axis' azimuths, at
several steps of a canopy's growth. And their geometry: a straightened tiller's stem reaches its straightening distance
(dredT) horizontally, then grows up. Without the option, AdelDyn's vertical axes (the comparison tests)."""
import numpy as np
import pytest

from openalea.adel.development import Planting
from openalea.adel.comparison.port_run import axis_names, organ_names, port_scene

PLANTS = [Planting(), Planting(x=12., azimuth=70., delay=40.)]


def get_inc_t(axis_name, haun_stage, inc_base, start_inc_t=1, inc_t_rate=30):
    """Adel.R's getIncT."""
    if axis_name == "MS":
        return inc_base
    if haun_stage > start_inc_t:
        return max(3, min(inc_base, inc_t_rate * (haun_stage - start_inc_t)))
    return 3


def axe_inclination(internode_visible, sheath_visible, inc_t, dred_t, epsilon=1e-6):
    """Adel.R's axe_inclination: (Einc, Ginc) per phytomer, from the visible lengths (Ev, Gv)."""
    number = len(internode_visible)
    e_inc, g_inc = [0.] * number, [0.] * number
    e_inc[0] = inc_t
    if dred_t > 0 and sum(internode_visible) + sum(sheath_visible) > epsilon:
        stem = []                                      # stemElements: (metamer, element, dl)
        for i in range(number):
            if internode_visible[i] > 0:
                stem.append((i, "en", internode_visible[i]))
            if sheath_visible[i] > 0:
                stem.append((i, "ga", sheath_visible[i]))
        alpha = (90 - inc_t) * np.pi / 180
        hc = np.cumsum([dl for _, _, dl in stem])
        dc = hc * np.cos(alpha)

        def set_inclination(element, value):
            metamer, kind, _ = stem[element]
            (e_inc if kind == "en" else g_inc)[metamer] = value

        if any(dc >= dred_t):
            nd = int(np.min(np.flatnonzero(dc >= dred_t)))        # 0-based
            if nd > 0:
                beta = np.arccos((dred_t - dc[nd - 1]) / (hc[nd] - hc[nd - 1]))
                set_inclination(nd, -(beta - alpha) / np.pi * 180)
                if nd < len(stem) - 1:
                    set_inclination(nd + 1, -(np.pi / 2 - beta) / np.pi * 180)
            else:
                beta = np.arccos(dred_t / hc[0])
                e_inc[0] = (np.pi / 2 - beta) / np.pi * 180
                if len(stem) > 1:
                    set_inclination(1, -e_inc[0])
    return e_inc, g_inc


@pytest.fixture(scope="module")
def grown():
    scene = port_scene("echap", plants=PLANTS, n_sect=2, parameters=dict(tiller_inclination=1))
    records = []
    for step in range(121):
        scene.run()
        if step % 20 == 0:
            ds = scene.populations[0].data_structure
            records.append({name: np.array(ds.get(name)) for name in (
                "inclination", "organ_azimuth", "axis_azimuth", "visible_length", "organ_rank", "haun_stage",
                "table_inclination", "straightening_distance", "insertion_azimuth", "first_azimuth")}
                           | dict(axes=axis_names(scene), organs=organ_names(ds),
                                  organ_axis=ds.map_to(np.arange(len(ds.entity_ids("Axis")), dtype=float), "Axis",
                                                       "Organ", "broadcast").astype(int)))
    return scene, records


def test_the_inclinations_of_adel_r(grown):
    _, records = grown
    compared = 0
    for record in records:
        for axis, name in enumerate(record["axes"]):
            mine = record["organ_axis"] == axis
            ranks = int(record["organ_rank"][mine].max())

            def by_rank(variable, organ):
                values = np.zeros(ranks)
                where = np.flatnonzero(mine & (record["organs"] == organ))
                values[record["organ_rank"][where] - 1] = record[variable][where]
                return values

            inc_t = get_inc_t(name, record["haun_stage"][axis], record["table_inclination"][axis])
            e_inc, g_inc = axe_inclination(by_rank("visible_length", "internode"), by_rank("visible_length", "sheath"),
                                           inc_t, record["straightening_distance"][axis])
            np.testing.assert_allclose(by_rank("inclination", "internode"), e_inc, atol=1e-9, err_msg=name)
            np.testing.assert_allclose(by_rank("inclination", "sheath"), g_inc, atol=1e-9, err_msg=name)
            assert record["axis_azimuth"][axis] == record["insertion_azimuth"][axis]
            assert by_rank("organ_azimuth", "internode")[0] == record["first_azimuth"][axis]
            compared += 1
    assert compared > 30


def test_the_straightened_tillers_reach_their_distance(grown):
    scene, records = grown
    ds = scene.populations[0].data_structure
    element_kind, _ = scene.populations[0].instance.architecture._kinds()
    node_axis = np.asarray(ds.owner("Axis"))
    x1, y1, x2, y2 = (np.asarray(ds.get(name)) for name in ("x1", "y1", "x2", "y2"))
    distance = np.asarray(ds.get("straightening_distance"))
    straightened = 0
    for axis, name in enumerate(axis_names(scene)):
        stem = np.flatnonzero((node_axis == axis) & (element_kind == 1))
        if name == "MS" or not (records[-1]["inclination"][records[-1]["organ_axis"] == axis] < 0).any():
            continue                                   # not straightened yet
        reach = np.hypot(x2[stem] - x1[stem[0]], y2[stem] - y1[stem[0]]).max()
        assert reach == pytest.approx(distance[axis], rel=0.1), name    # the rolled blades stack on it too
        straightened += 1
    assert straightened >= 6


def test_vertical_without_the_option():
    scene = port_scene("echap", plants=PLANTS, n_sect=2)
    for _ in range(80):
        scene.run()
    ds = scene.populations[0].data_structure
    assert not np.any(np.asarray(ds.get("inclination")))
    assert not np.any(np.asarray(ds.get("axis_azimuth")))
