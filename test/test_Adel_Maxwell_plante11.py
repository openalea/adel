"""Test stability of Adel for simulating a reference plant (Maxwell, plant 11, pareameterised by M. Abichou/B.Andrieu, EGC Grignon"""
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



import pathlib
from functools import reduce

from openalea.adel.AdelR import devCsv, setAdel, RunAdel, genGeoLeaf, genGeoAxe, \
    csvAsDict


def test_organ_length():
    dir = str(pathlib.Path(__file__).parent.resolve() / "data" / "test_Adel_Maxwell_plante11/Maxwell_")
    sufix = "_plante11.csv"
    axeTpath = dir + "axeT" + sufix
    dimTpath = dir + "dimT" + sufix
    phenTpath = dir + "phenT" + sufix
    earTpath = dir + "earT.csv"
    ssi2senpath = dir + "ssi2sen.csv"
    devT = devCsv(axeTpath, dimTpath, phenTpath, earTpath, ssi2senpath)
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT, geoLeaf, geoAxe, 1, seed=1)
    cantables = [RunAdel(x, pars) for x in range(0, 2300, 100)]
    expected = csvAsDict(dir + "reference_simulation.csv")
    for k in (
        "TT",
        "plant",
        "numphy",
        "Ll",
        "Gl",
        "El",
        "Lv",
        "Gv",
        "Ev",
        "Esen",
        "Lsen",
        "Gsen",
    ):
        sim = reduce(lambda x, y: x + y, (t[k].tolist() for t in cantables))
        exp = expected[k]
        # np.testing.assert_allclose(sim,exp)
