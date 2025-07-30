"""adel dataflow Tests"""

__license__ = "Cecill-C"
__revision__ = " $Id: $"

from openalea.adel.sensitivity.sensitivity import *


def test1():
    repet = 5
    factors = ["a", "b"]
    binf = (0, 0)
    bsup = (1, 1)

    m, pdict = Morris(repet, factors, binf, bsup)
    return m, pdict
