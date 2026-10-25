"""ShootArchitecture on a stand of four plants in one MPG against the legacy's stand (AdelDyn.build_stand), plant by
plant: each at its position and azimuth, its development delayed (the plants differ). Every drawn element's anchor
point and triangles, and the elements' lengths and areas.

The legacy places each plant twice (the turtle, then update_geometry's transform_geom: its position again, its
azimuth read in radians); its triangles are compared with that second placement undone, and its T1 offset taken off.
"""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import Planting, legacy_run
from openalea.adel.comparison.port_run import port_run

from test_geometry import T1_OFFSET, legacy_drawn

STEPS, EVERY = 120, 30
PLANTS = [Planting(0., 0., 0., 0.), Planting(15., 0., 40., 30.), Planting(0., 12., 130., 60.),
          Planting(15., 12., 250., 15.)]


@pytest.fixture(scope="module")
def runs():
    return (legacy_run(n_steps=STEPS, leaves="echap", record_every=EVERY, n_sect=2, plants=PLANTS),
            port_run(n_steps=STEPS, leaves="echap", record_every=EVERY, n_sect=2, plants=PLANTS))


def turtle_triangles(triangles, plant, axis):
    """The legacy's triangles (cm) without update_geometry's second placement and T1 offset."""
    angle = -plant.azimuth                             # the legacy gives the degrees to an angle in radians
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.], [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    out = (triangles - [plant.x, plant.y, 0.]) @ rotation.T
    if axis.startswith("T1"):
        out[..., 2] -= T1_OFFSET
    return out


def test_every_plant_as_the_legacy(runs):
    legacy, port = runs
    compared = np.zeros(len(PLANTS), dtype=int)
    for legacy_record, port_records in zip(legacy.steps, port.steps):
        for k, plant in enumerate(PLANTS):
            for key, label, element in legacy_drawn(legacy_record, plant=k + 1):
                values = port_records[k][key]["elements"][label]
                anchor = np.array(list(legacy_record.props["anchor_point"][element]), dtype=float) * 1e-2
                assert [values["x1"], values["y1"], values["z1"]] == pytest.approx(anchor, abs=1e-12), (k, key, label)
                assert values["length"] == pytest.approx(legacy_record.props["length"][element] * 1e-2, abs=1e-12)
                expected = turtle_triangles(legacy_record.triangles[element], plant, key[0]) * 1e-2
                assert values["triangles"].shape == expected.shape, (k, key, label)
                assert values["triangles"] == pytest.approx(expected, abs=1e-10), (k, key, label)
                compared[k] += 1
    assert (compared > 50).all(), compared


def test_the_plants_differ(runs):
    _, port = runs
    areas = [sum(organ["organ_area"] for organ in record.values()) for record in port.steps[1]]
    assert len(set(np.round(areas, 12))) == len(PLANTS)             # delayed plants are smaller
