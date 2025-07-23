import openalea.adel.mtg as mtg
import openalea.adel.fitting as fitting
import openalea.adel.parameterisation.parameterisation as adelp
from openalea.adel.symbol import build_symbols
from openalea.adel.data_samples import leaves_db as ldb
from openalea.mtg.traversal import pre_order2
from openalea.plantgl.all import Viewer


g = None


def adel_mtg():
    global g
    if g is not None:
        return g

    pars = adelp.simpleMais_param(
        total_area=5000,
        total_height=200,
        pseudostem_height=20,
        nb_phy=15,
        nb_young_phy=6,
        lad_skew=5,
        lad_rmax=0.7,
        pseudostem_dist=1.4,
        stem_dist=1.0,
        phyllotactic_angle=180,
        basal_insertion=40,
        top_insertion=30,
        diam_base=2,
        diam_top=1,
        leaf_width_ratio=0.1,
        shape_factors=0.75,
    )
    canT = adelp.simpleMais2dict(pars)[0]
    g = mtg.mtg_factory(canT, sectors=1)

    return g


def leaves_db():
    leaves = ldb()
    leaves = fitting.fit_leaves(leaves, 9)

    db = leaves[0]
    functions = build_symbols(db)
    return functions


def test_adel_mtg():
    g = adel_mtg()
    # print g
    symbols = leaves_db()

    assert len(g) == len(g.sub_mtg(g.root))

    print(list(symbols.keys()))
    g = mtg.mtg_turtle(g, symbols)
    assert g is not None


def test_dynamic_mtg():
    """Add to each metamer a thermal time (e.g. 10).
    Length = length * (thermal_time - start_thermaltime)/(end-start)
    if the thermal time is lower that the start thermal time,
    do not proceed more.
    """

    g = adel_mtg()

    v = next(g.component_roots_at_scale_iter(g.root, scale=3))
    tt = 0
    dtt = 10.0
    for metamer in pre_order2(g, v):
        nm = g.node(metamer)
        for node in nm.components():
            node.start_tt = tt
            node.end_tt = tt + dtt
        tt += dtt

    assert g is not None and tt != 0


def display():
    g, max_time = test_dynamic_mtg()
    symbols = leaves_db()

    for time in range(1, 150):
        g.properties()["geometry"] = {}
        g = mtg.mtg_turtle_time(g, symbols, time)
        scene = g.to_plantgl()
        Viewer.display(scene)
