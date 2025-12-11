import pytest

try:
    # test if R's stats package loads correctly
    import rpy2.robjects as ro
    ro.r("na.omit")  # will error on broken R
    R_OK = True
except Exception:
    R_OK = False

pytest_r_skip = pytest.mark.skipif(
    not R_OK,
    reason="Skipping R/rpy2 tests on this platform due to broken R/rpy2 installation"
)
