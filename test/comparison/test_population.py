"""ShootArchitecture on a stand of four plants in one MPG against the legacy's stand (AdelDyn.build_stand), plant by
plant: each at its position and azimuth, its development delayed (the plants differ). Every drawn element's anchor
point and triangles, and the elements' lengths and areas.

The legacy places each plant twice (the turtle, then update_geometry's transform_geom: its position again, its
azimuth read in radians); its triangles are compared with that second placement undone, and its T1 offset taken off.
The port starts T1 at the coleoptile's top: the legacy's T1 is moved there (port_run.coleoptile_top).
"""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import StandPlant, legacy_run, turtle_placement
from openalea.adel.comparison.port_run import coleoptile_top, port_run

from test_geometry import legacy_drawn

STEPS, EVERY = 120, 30
PLANTS = [StandPlant(0., 0., 0., 0.), StandPlant(15., 0., 40., 30.), StandPlant(0., 12., 130., 60.),
          StandPlant(15., 12., 250., 15.)]


@pytest.fixture(scope="module")
def runs():
    return (legacy_run(n_steps=STEPS, leaves="echap", record_every=EVERY, n_sect=2, plants=PLANTS),
            port_run(n_steps=STEPS, leaves="echap", record_every=EVERY, n_sect=2, plants=PLANTS))


def test_every_plant_as_the_legacy(runs):
    legacy, port = runs
    compared = np.zeros(len(PLANTS), dtype=int)
    for legacy_record, port_records in zip(legacy.steps, port.steps):
        for k, plant in enumerate(PLANTS):
            for key, label, element in legacy_drawn(legacy_record, plant=k + 1):
                values = port_records[k][key]["elements"][label]
                anchor = np.array(list(legacy_record.props["anchor_point"][element]), dtype=float) * 1e-2
                top = coleoptile_top(port_records[k])          # the port starts T1 at the coleoptile's top
                if key[0].startswith("T1"):
                    anchor += top
                assert [values["x1"], values["y1"], values["z1"]] == pytest.approx(anchor, abs=1e-12), (k, key, label)
                assert values["length"] == pytest.approx(legacy_record.props["length"][element] * 1e-2, abs=1e-12)
                expected = turtle_placement(legacy_record.triangles[element], plant, key[0], top * 1e2) * 1e-2
                assert values["triangles"].shape == expected.shape, (k, key, label)
                assert values["triangles"] == pytest.approx(expected, abs=1e-10), (k, key, label)
                compared[k] += 1
    assert (compared > 50).all(), compared


def test_the_plants_differ(runs):
    _, port = runs
    areas = [sum(organ["organ_area"] for organ in record.values()) for record in port.steps[1]]
    assert len(set(np.round(areas, 12))) == len(PLANTS)             # delayed plants are smaller
