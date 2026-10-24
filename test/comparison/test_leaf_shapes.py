"""LeafShapes' fitted midribs and age classes against the legacy's (geometric_elements.Leaves: fitting.fit_leaves,
get_age_index), for Adel's samples and echap's database."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import leaf_database
from openalea.adel.comparison.port_run import leaf_shapes


@pytest.mark.parametrize("name", ["samples", "echap"])
def test_the_same_fitted_midribs(name):
    legacy, port = leaf_database(name), leaf_shapes(name)
    legacy = legacy[0] if isinstance(legacy, dict) else legacy
    compared = 0
    for shape_class, key in enumerate(sorted(legacy.leaves, key=int), start=1):
        for seed, shape in enumerate(legacy.leaves[key], start=1):
            for age_class, fitted in (shape.items() if isinstance(shape, dict) else [(0, shape)]):
                index = port.shape_index([shape_class], [seed], [age_class])[0]
                assert index >= 0
                for legacy_values, port_values in zip(fitted, port.fitted[index]):
                    assert port_values == pytest.approx(np.asarray(legacy_values), abs=1e-12)
                compared += 1
    assert compared == len(port.fitted) > 0


def test_the_same_age_classes():
    legacy, port = leaf_database("echap"), leaf_shapes("echap")
    ages = np.linspace(-3., 20., 231)
    assert port.age_class(ages).tolist() == [legacy.get_age_index(age) for age in ages]
    assert leaf_shapes("samples").age_class(ages).tolist() == [0] * ages.size
