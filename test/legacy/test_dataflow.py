"""adel dataflow Tests"""

__license__ = "Cecill-C"
__revision__ = " $Id: $"

from openalea.core.alea import run
from openalea.core.pkgmanager import PackageManager
import warnings
from random import random, randint

""" A unique PackageManager is created for all test of dataflow """
pm = PackageManager()
pm.init(verbose=False)


def test_adelr1():
    """Test AdelR MonoRun"""
    # check if openalea components are there
    if "#openalea.python" in [k for k, v in pm.items()]:
        res = run(
            ("openalea.adel.tutorials", "AdelR MonoRun"), inputs={}, pm=pm, vtx_id=13
        )
    else:
        warnings.warn(
            "OpenAlea components not installed, AdelR Monorun dataflow could not be test"
        )


# to be repaired
# def test_adelr2():
#     """ Test AdelR MonoRun """
#     res = run(('openalea.adel.tutorials', 'AdelR MonoRun'),
#         inputs={}, pm=pm, vtx_id=39)
#
#
# def test_arvalis():
#     """ Test Leaf db Explorer """
#     res = run(('openalea.adel.tutorials', 'Leaf db Explorer'),
#         inputs={}, pm=pm)
#     assert res[0] == []
