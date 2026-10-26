"""ShootArchitecture against the legacy (AdelDyn) on plants the harness' sample plant does not have: a tiller of a
tiller, T1.1 (on T1's phytomer 1, its first leaf after T4's: the axes table in its order of emergence, as the port
reads it), and blades without a shape (the phytomer table's Lindex or Lseed -999: the legacy's shape key None). The
same comparisons as test_structure, test_elements and test_geometry: the organs, the elements' lengths and areas, the
anchor points and the triangles.

The legacy's turtle (AD-37) continues T1 after T1.1 on the main stem's top saved, and starts T1.1 on the tillers'
base; the port reproduces it by default. With own_axis_top=1 (the correction), a tiller of a tiller starts on its
bearer's phytomer and its bearer continues on its own top: T1 is placed as without T1.1, and a plant without tillers
of tillers is placed as by Adel's turtle."""
import os

import numpy as np
import pandas as pd
import pytest

from openalea.adel.comparison.legacy_run import T1_OFFSET, legacy_run
from openalea.adel.comparison.port_run import coleoptile_top, port_run
from openalea.adel.tables import DATA, with_tillers

import test_elements as elements
import test_geometry as geometry
from test_structure import legacy_structure

STEPS, EVERY = 120, 20
INCLINATIONS = {("MS", 3): 10., ("T1", 1): 30., ("T1", 2): 15., ("T1.1", 1): 35., ("T1.1", 2): 12., ("T2", 1): 25.}
SCENARIOS = {"vertical": dict(plant_azimuth=0., inclinations=None),
             "inclined, azimuth 0": dict(plant_azimuth=0., inclinations=INCLINATIONS),
             "inclined": dict(plant_azimuth=37., inclinations=INCLINATIONS)}
SHAPELESS = (("MS", 4, "blade"), ("T2", 2, "blade"))


def tables(tiller_of_tiller=True):
    """(axes table, phytomer table): Adel's sample plant, with T1.1 (T3's row, 5 leaves, its first at 420 °Cd), the
    main stem's blade 4 without a shape class and T2's blade 2 without a shape seed."""
    axes = pd.read_csv(os.path.join(DATA, "axeTable.csv"))
    if tiller_of_tiller:
        tiller = axes.loc[axes["axe"] == "T3"].assign(axe="T1.1", emf1=420., nf=5., nf_end=5.)
        axes = pd.concat([axes, tiller], ignore_index=True).sort_values("emf1", kind="stable").reset_index(drop=True)
    phytomers = with_tillers(axes, pd.read_csv(os.path.join(DATA, "phytoT.csv")))
    phytomers.loc[(phytomers["axe"] == "MS") & (phytomers["n"] == 4), "Lindex"] = -999.
    phytomers.loc[(phytomers["axe"] == "T2") & (phytomers["n"] == 2), "Lseed"] = -999.
    return axes, phytomers


def run_options(name):
    return dict(n_steps=STEPS, record_every=EVERY, leaves="samples", n_sect=3, **SCENARIOS[name])


@pytest.fixture(scope="module", params=list(SCENARIOS))
def runs(request):
    options = run_options(request.param)
    return request.param, (legacy_run(tables=tables(), **options), port_run(tables=tables(), **options))


def test_the_same_organs(runs):
    _, (legacy, port) = runs
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        assert set(port_record) == set(legacy_structure(legacy_record, None))
    assert any(key[0] == "T1.1" for key in port.steps[-1])


def test_the_same_elements(runs):
    _, (legacy, port) = runs
    elements.test_the_same_elements((("samples", 3), legacy, port))


def test_the_same_anchor_points(runs):
    _, (legacy, port) = runs
    geometry.test_the_same_anchor_points((None, (legacy, port)))


def test_the_same_triangles(runs):
    name, (legacy, port) = runs
    if SCENARIOS[name]["plant_azimuth"] != 0.:
        pytest.skip("the legacy places the plant twice: its triangles are compared at the azimuth 0")
    compared = 0
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        for key, label, element in geometry.legacy_drawn(legacy_record):
            expected = legacy_record.triangles[element] * 1e-2
            if key[0].startswith("T1"):
                expected = expected - T1_OFFSET * 1e-2 * np.array([0., 0., 1.]) + coleoptile_top(port_record)
            triangles = port_record[key]["elements"][label]["triangles"]
            assert triangles.shape == expected.shape and triangles == pytest.approx(expected, abs=1e-10), (key, label)
            compared += 1
    assert compared > 200


def test_the_blades_without_a_shape(runs):
    """Their sectors have the legacy's lengths and flat parts, no area, no width and no triangles."""
    _, (legacy, port) = runs
    legacy_record, port_record = legacy.steps[-1], port.steps[-1]
    drawn = {key for key, _, _ in geometry.legacy_drawn(legacy_record)}
    for key in SHAPELESS:
        organ = port_record[key]
        assert -1 in (organ["shape_class"], organ["shape_seed"]) and organ["organ_area"] == 0. and key not in drawn
        leaf = [values for label, values in organ["elements"].items() if label.startswith("LeafElement")]
        assert sum(values["length"] for values in leaf) > 0.
        assert all(values["area"] == 0. and values["width"] == 0. and len(values["triangles"]) == 0
                   for values in leaf)


def positions(run):
    """{(axis, rank, organ, element): (x1, y1, z1)} of a port run's last step."""
    return {key + (label,): np.array([values["x1"], values["y1"], values["z1"]])
            for key, organ in run.steps[-1].items() for label, values in organ["elements"].items()}


@pytest.mark.parametrize("name", ["vertical", "inclined"])
def test_own_axis_top(name):
    """With own_axis_top=1, T1.1 starts on its bearer's element (T1's first internode's first) and T1 is placed as
    without T1.1; without it (Adel's turtle), T1 after T1.1 is not."""
    corrected = positions(port_run(tables=tables(), parameters=dict(own_axis_top=1), **run_options(name)))
    without = positions(port_run(tables=tables(tiller_of_tiller=False), parameters=dict(own_axis_top=1),
                                 **run_options(name)))
    for key, position in without.items():
        assert corrected[key] == pytest.approx(position, abs=1e-12), key
    assert corrected[("T1.1", 1, "internode", "HiddenElement")] == pytest.approx(
        corrected[("T1", 1, "internode", "HiddenElement")], abs=1e-12)
    if name == "inclined":
        legacy = positions(port_run(tables=tables(), **run_options(name)))
        assert not all(np.allclose(legacy[key], position, atol=1e-9) for key, position in without.items()
                       if key[0] == "T1")


@pytest.mark.parametrize("name, tiller_inclination", [("vertical", 0), ("inclined", 0), ("vertical", 1)])
def test_own_axis_top_changes_no_ordinary_plant(name, tiller_inclination):
    """On the harness' plant (no tiller of a tiller), the correction places every element as Adel's turtle, with
    the given inclinations or static Adel's (tiller_inclination=1)."""
    options = run_options(name)
    corrected = positions(port_run(parameters=dict(own_axis_top=1, tiller_inclination=tiller_inclination),
                                   **options))
    for key, position in positions(port_run(parameters=dict(tiller_inclination=tiller_inclination),
                                            **options)).items():
        assert corrected[key] == pytest.approx(position, abs=1e-12), key
