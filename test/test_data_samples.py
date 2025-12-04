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

from openalea.adel.data_samples import *


def test_leaf_db():
    db = leaves_db()
    sr = srdb()
    xy = xydb()
    leaves = wheat_leaf_db()
    assert db is not None
    assert sr is not None
    assert xy is not None
    assert leaves is not None


def test_devT():
    dev = devT()
    assert dev is not None


def test_adel():
    two_metamer = adel_two_metamers()
    two_metamer_2 = adel_two_metamers(leaf_sectors=2)
    one_leaf = adel_one_leaf()
    one_sect = adel_one_leaf_element()
    assert two_metamer is not None
    assert two_metamer_2 is not None
    assert one_leaf is not None
    assert one_sect is not None


def test_stand():
    g, domain_area, domain, convUnit = adel_two_metamers_stand()
    assert g is not None
    assert domain_area is not None
    assert domain is not None
    assert convUnit is not None
