"""ShootArchitecture's elements against the legacy's (AdelDyn.update_geometry's update_organ_elements), on the
harness' scenario: every element's lengths (green, senesced, rolled, hidden), its flat part's position and width
(1e-9), and its areas and the organs' visible areas (the legacy in cm, the port in m). The areas differ: the port
integrates the blades' width profiles exactly (linear between their points), the legacy by Simpson's rule on those
points, least exact on small parts (a tip just out): every area within 8 %, their median within 2e-3.

The widths are the mean of the width profile at a part's bounds and its points strictly inside. With echap's
profiles (points every 0.05) and 4 elements per blade, the elements' bounds 0.25, 0.5, 0.75 fall on points: the
legacy computes them in cm, and some round one ulp below (0.7499999999999999), the point then inside the element
above; the port counts a point within 1e-9 of a bound as the bound (the definition's exact result). Those widths
differ (a sixth of them, by less than 10 %); every other is the legacy's to 1e-9."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import legacy_run
from openalea.adel.comparison.port_run import port_run

STEPS, EVERY = 120, 20
LENGTHS = ("length", "green_length", "senesced_length")
AREAS = ("area", "green_area", "senesced_area")


def legacy_elements(record):
    """{(axis, rank, organ): ({element label: {variable: value}}, organ area)} in m and m2."""
    out = {}
    for organ in record.vertices:
        if record.scale[organ] != 4 or record.label[organ] not in ("internode", "sheath", "blade"):
            continue
        metamer = record.complex[organ]
        key = (record.label[record.complex[metamer]], int(record.label[metamer][len("metamer"):]),
               record.label[organ])
        elements = {}
        for e in record.vertices:
            if record.scale[e] == 5 and record.complex[e] == organ and record.label[e] not in ("baseElement",
                                                                                               "topElement"):
                props = {name: values[e] for name, values in record.props.items() if e in values}
                elements[record.label[e]] = props
        out[key] = (elements, (record.props["area"].get(organ) or 0.) * 1e-4)
    return out


BOUNDS_ON_POINTS = {("echap", 4)}                     # the elements' bounds on the width profiles' points


@pytest.fixture(scope="module", params=[("samples", 1), ("echap", 1), ("echap", 4), ("samples", 3)])
def runs(request):
    leaves, n_sect = request.param
    return (request.param, legacy_run(n_steps=STEPS, leaves=leaves, record_every=EVERY, n_sect=n_sect),
            port_run(n_steps=STEPS, leaves=leaves, record_every=EVERY, n_sect=n_sect))


def test_the_same_elements(runs):
    run_name, legacy, port = runs
    differences, width_differences = [], []
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        for key, (elements, organ_area) in legacy_elements(legacy_record).items():
            port_elements = port_record[key]["elements"]
            for label, legacy_values in elements.items():
                values = port_elements[label]
                for name in LENGTHS:
                    if legacy_values.get(name) is not None:
                        assert values[name] == pytest.approx(legacy_values[name] * 1e-2, rel=1e-9, abs=1e-12), (
                            key, label, name)
                for name in AREAS:
                    expected = legacy_values.get(name)
                    if expected:
                        differences.append(abs(values[name] - expected * 1e-4) / (expected * 1e-4))
                    elif expected is not None:
                        assert values[name] == pytest.approx(0., abs=1e-12), (key, label, name)
                if label.startswith("LeafElement"):
                    assert values["rolled_part"] == pytest.approx(legacy_values["lrolled"] * 1e-2, abs=1e-12)
                    for port_name, legacy_name in (("flat_base", "srb"), ("flat_top", "srt")):
                        expected = legacy_values[legacy_name]
                        if expected is None:
                            assert np.isnan(values[port_name]), (key, label, port_name)
                        else:
                            assert values[port_name] == pytest.approx(expected, rel=1e-9, abs=1e-12)
                    if run_name in BOUNDS_ON_POINTS and legacy_values["width"]:
                        expected = legacy_values["width"] * 1e-2
                        width_differences.append(abs(values["width"] - expected) / expected)
                    else:
                        assert values["width"] == pytest.approx(legacy_values["width"] * 1e-2, rel=1e-9,
                                                                abs=1e-12)
                assert bool(values["is_green"]) == bool(legacy_values["is_green"]), (key, label)
            if organ_area > 0:
                differences.append(abs(port_record[key]["organ_area"] - organ_area) / organ_area)
    if width_differences:
        width_differences = np.array(width_differences)
        assert np.mean(width_differences > 1e-9) < 0.2 and width_differences.max() < 0.1
    differences = np.array(differences)
    assert differences.size and differences.max() < 0.08 and np.median(differences) < 2e-3, (
        differences.max(), np.median(differences))
