"""Test use cases of adel/fspm wheat coupling"""
import sys
import pytest
import rpy2

try:
    # test if R's stats package loads correctly
    import rpy2.robjects as ro
    ro.r("na.omit")  # will error on broken R
    r_ok = True
except Exception:
    r_ok = False

pytestmark = pytest.mark.skipif(
    sys.platform == "win32" or not r_ok,
    reason="Skipping R tests on Windows CI due to broken R installation"
)

from openalea.adel.adelwheat_dynamic import AdelWheatDyn


def test_add_metamer():
    # create a plant with 3 metamers
    adel = AdelWheatDyn(seed=1234)
    g = adel.setup_canopy()
    labels = g.property("label")
    metamers = [vid for vid in labels if labels[vid].startswith("metamer")]
    assert len(metamers) == 3

    # add empty new metamer
    vid = adel.add_metamer(g, 1, "MS")
    new_metamer = g.node(vid)
    metamers = [vid for vid in labels if labels[vid].startswith("metamer")]
    assert len(metamers) == 4
    internode, sheath, blade = new_metamer.components()
    assert blade.length == 0
    elts = [c.label for c in blade.components()]
    assert elts == ["baseElement", "topElement"]

    # grow leaf and check components
    blade.length = 6
    blade.visible_length = 3
    adel.update_geometry(g)
    assert blade.area > 0
    elts = [c.label for c in blade.components()]
    assert len(elts) > 2
