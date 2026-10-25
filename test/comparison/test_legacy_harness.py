"""The legacy harness (openalea.adel.comparison.legacy_run): deterministic, its records complete, the visible
elements meshed."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import legacy_run


@pytest.fixture(scope="module", params=["samples", "echap"])
def runs(request):
    return [legacy_run(n_steps=80, leaves=request.param, record_every=40) for _ in range(2)]


def test_deterministic(runs):
    first, second = runs
    for a, b in zip(first.steps, second.steps):
        assert a.vertices == b.vertices and a.label == b.label and a.parent == b.parent
        for name, values in a.props.items():
            for v, value in values.items():
                other = b.props[name][v]
                if isinstance(value, float):
                    assert value == other or (np.isnan(value) and np.isnan(other)), name
        assert set(a.triangles) == set(b.triangles)
        for v in a.triangles:
            np.testing.assert_array_equal(a.triangles[v], b.triangles[v])


def test_the_main_stem_grows(runs):
    steps = runs[0].steps
    metamers = [len([v for v in s.vertices if s.scale[v] == 3 and s.label[v] != "metamer0"]) for s in steps]
    assert metamers == sorted(metamers) and metamers[-1] > metamers[0]


def test_visible_elements_are_meshed(runs):
    last = runs[0].steps[-1]
    length = last.props["length"]
    visible = [v for v in last.vertices
               if last.scale[v] == 5 and last.label[v].startswith(("StemElement", "LeafElement"))
               and length.get(v, 0) > 0]
    assert visible
    meshed = [v for v in visible if v in last.triangles and last.triangles[v].shape[0] > 0]
    assert len(meshed) >= 0.9 * len(visible)
    for v in meshed:
        assert np.isfinite(last.triangles[v]).all()
