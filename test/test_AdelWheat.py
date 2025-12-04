import pytest

try:
    # test if R's stats package loads correctly
    import rpy2.robjects as ro
    ro.r("na.omit")  # will error on broken R
    r_ok = True
except Exception:
    r_ok = False

pytestmark = pytest.mark.skipif(
    not r_ok,
    reason="Skipping R tests on Windows CI due to broken R installation"
)

from openalea.adel.astk_interface import AdelWheat
from openalea.astk.Weather import sample_weather




def test_adelwheat():
    seq, weather = sample_weather()
    seq = seq.tz_localize(tz='Europe/Paris')
    wdata = weather.get_weather(seq)

    adel = AdelWheat(nsect=2)

    g = adel.setup_canopy(100)
    adel.grow(g, wdata)
