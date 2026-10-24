"""ShootArchitecture's placement against the legacy's turtle (mtg_interpreter.AdelVisitor), on the harness' scenario:
every drawn element's anchor point (the legacy's anchor_point), and the stem elements' triangles (the legacy in cm,
the port in m).

Two scenarios: the legacy's (no inclination, the plant's azimuth 0), and a plant azimuth with inclined internodes
(the main stem's 3rd, the tillers' first: the turtle's rules for an axis' first inclined element and for the others).

Known difference: AdelDyn.update_geometry translates the elements of the axes named T1… by 0.0285 cm up (a fix of
the first tiller's position) and places the plant a second time (its position and azimuth, the azimuth in radians);
the port places them once, by the turtle. The triangles are compared with the plant's azimuth 0, the legacy's T1
offset taken off."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import legacy_run
from openalea.adel.comparison.port_run import port_run

from test_elements import legacy_elements

STEPS, EVERY = 120, 20
T1_OFFSET = 0.0285                                     # cm, AdelDyn.update_geometry's
INCLINATIONS = {("MS", 3): 10., ("T1", 1): 30., ("T2", 1): 25., ("T1", 2): 15., ("T3", 1): 20.}
SCENARIOS = {"vertical": dict(plant_azimuth=0., inclinations=None),
             "inclined": dict(plant_azimuth=37., inclinations=INCLINATIONS),
             "inclined, azimuth 0": dict(plant_azimuth=0., inclinations=INCLINATIONS)}


def legacy_drawn(record):
    """[(step's organ key, element label, element vid)] of the record's drawn elements (with an anchor point)."""
    anchors = record.props.get("anchor_point", {})
    organs = {}
    for organ in record.vertices:
        if record.scale[organ] == 4 and record.label[organ] in ("internode", "sheath", "blade"):
            metamer = record.complex[organ]
            organs[organ] = (record.label[record.complex[metamer]], int(record.label[metamer][len("metamer"):]),
                             record.label[organ])
    return [(organs[record.complex[e]], record.label[e], e) for e in record.vertices
            if record.scale[e] == 5 and e in anchors and record.complex[e] in organs]


@pytest.fixture(scope="module", params=list(SCENARIOS))
def runs(request):
    scenario = SCENARIOS[request.param]
    return request.param, (legacy_run(n_steps=STEPS, leaves="echap", record_every=EVERY, **scenario),
                           port_run(n_steps=STEPS, leaves="echap", record_every=EVERY, **scenario))


def test_the_same_anchor_points(runs):
    _, (legacy, port) = runs
    compared = 0
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        for key, label, element in legacy_drawn(legacy_record):
            values = port_record[key]["elements"][label]
            expected = np.array(list(legacy_record.props["anchor_point"][element]), dtype=float) * 1e-2
            assert [values["x1"], values["y1"], values["z1"]] == pytest.approx(expected, abs=1e-12), (key, label)
            compared += 1
    assert compared > 200


def test_the_same_stem_triangles(runs):
    name, (legacy, port) = runs
    if SCENARIOS[name]["plant_azimuth"] != 0.:
        pytest.skip("the legacy places the plant twice: its triangles are compared at the azimuth 0")
    compared = 0
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        for key, label, element in legacy_drawn(legacy_record):
            if label != "StemElement":
                continue
            expected = legacy_record.triangles[element].copy()
            if key[0].startswith("T1"):
                expected[..., 2] -= T1_OFFSET
            assert port_record[key]["elements"][label]["triangles"] == pytest.approx(expected * 1e-2, abs=1e-12), (
                key, label)
            compared += 1
    assert compared > 100
