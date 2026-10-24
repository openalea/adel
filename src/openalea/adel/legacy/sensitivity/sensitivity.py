#
# Utilitaires pour l'analyse de sensibilite
#
import numpy
import rpy2.robjects as robj

r = robj.r

r.require("sensitivity")


def Morris(repeat, factors, binf, bsup):
    """Simplified import of R'Morris function"""

    factors = robj.StrVector(factors)
    binf = numpy.array(binf)
    bsup = numpy.array(bsup)
    d = r.list("oat", 5, 3)
    d = r["names<-"](d, ["type", "levels", "grid.jump"])
    m = r.morris(factors=factors, r=repeat, design=d, binf=binf, bsup=bsup)
    # param=r['data.frame'](m.rx["X"])
    param = r["data.frame"](r["$"](m, "X"))
    pdict = dict((k, numpy.array((r["$"](param, k))).tolist()) for k in param.names)
    # pdict = dict((str(k), list(v)) for k,v in param.iteritems())
    return m, pdict


def Morris_IS(morrisPlan, y):
    """Compute morris sensitivity indicators for a Morris plan and the corresponding model responses"""
    mod = r.tell(morrisPlan, numpy.array(y))
    ee = r["$"](mod, "ee")
    mu = r.apply(ee, 2, r["mean"])
    mustar = r.apply(ee, 2, r("function(x) mean(abs(x))"))
    sigma = r.apply(ee, 2, r.sd)

    return mod, numpy.array(mustar).tolist(), numpy.array(sigma).tolist()


def plotSens(mod):
    """moris plot"""
    return r.plot(mod)
