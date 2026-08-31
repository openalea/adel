"""Test stability of Adel for simulating a reference plant (Maxwell, plant 11, pareameterised by M. Abichou/B.Andrieu, EGC Grignon"""
from pathlib import Path
from functools import reduce
import numpy as np

from openalea.adel.AdelR import devCsv, csvAsDict
from openalea.adel.astk_interface import AdelWheat

from .conftest import pytest_r_skip


@pytest_r_skip
def test_organ_length():
    datadir = Path(__file__).parent if '__file__' in globals() else Path(".")
    dir = str(datadir.resolve() / "data" / "test_Adel_Maxwell_plante11/Maxwell_")
    sufix = "_plante11.csv"
    axeTpath = dir + "axeT" + sufix
    dimTpath = dir + "dimT" + sufix
    phenTpath = dir + "phenT" + sufix
    earTpath = dir + "earT.csv"
    ssi2senpath = dir + "ssi2sen.csv"
    devT = devCsv(axeTpath, dimTpath, phenTpath, earTpath, ssi2senpath)
    adel = AdelWheat(devT=devT, seed=1)
    cantables = [adel.run_adel(x) for x in range(0, 2300, 100)]
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
        np.testing.assert_allclose(sim,exp)
