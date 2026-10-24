from rpy import *


def Ecalc(eabs, area, table):
    """ """

    # write the node code here.

    tp = {"ind": list(range(len(eabs))), "eabs": eabs, "area": area}
    table.update(tp)

    r.load(
        "C:/Documents and Settings/Jessica/Mes documents/openalea_pkg/Eabscalc.RData"
    )
    res = r.Eabscalc(table)

    # return outputs
    return (res["area sum"], res["eabs moy"])
