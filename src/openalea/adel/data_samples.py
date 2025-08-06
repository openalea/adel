"""Provides sample data for testing Adel without reference to file paths, nor needs to call package manager"""

import os
import pandas
from openalea.adel.geometric_elements import Leaves, xydb_from_csv, srdb_from_csv
import openalea.adel.json_numpy as json_np

from openalea.adel.newmtg import mtg_factory, adel_metamer
from openalea.adel.mtg_interpreter import mtg_interpreter
from openalea.adel.stand.stand import regular


datadir = os.path.dirname(__file__)


def wheat_dimension_profiles():
    """Dimension profiles compiled using maxwell REF + Echap data (cf Dimension_profile.R in echap/src/data/architecturaldata)"""
    fn = datadir + "/data/Wheat_dimension_profiles.csv"
    df = pandas.read_csv(fn)
    return df


def leaves():
    return {0: Leaves()}


def leaves_db():
    import openalea.adel.fitting as fitting

    fn = datadir + "/data/simpleleavesdb.json"
    with open(fn) as f:
        leaves = json_np.load(f)
    leaves, discard = fitting.fit_leaves(leaves, 9)
    return leaves


def _devT():
    from openalea.adel.AdelR import _devCsv

    axeT = datadir + "/data/axeTCa0N.csv"
    dimT = datadir + "/data/dimTCa0N.csv"
    phenT = datadir + "/data/phenTCa0N.csv"
    earT = datadir + "/data/earTCa0N.csv"
    ssisenT = datadir + "/data/ssi2sen.csv"

    return _devCsv(axeT, dimT, phenT, earT, ssisenT)


def devT():
    from openalea.adel.AdelR import devCsv

    axeT = datadir + "/data/axeTCa0N.csv"
    dimT = datadir + "/data/dimTCa0N.csv"
    phenT = datadir + "/data/phenTCa0N.csv"
    earT = datadir + "/data/earTCa0N.csv"
    ssisenT = datadir + "/data/ssi2sen.csv"

    return devCsv(axeT, dimT, phenT, earT, ssisenT)


def srdb(nlevel=None):
    rdata = datadir + "/data/SRSo.csv"
    db = srdb_from_csv(rdata)
    if nlevel is None:
        return db
    else:
        return {str(k + 1): db["1"] for k in range(nlevel)}


def xydb():
    data = datadir + "/data/So99.csv"
    return xydb_from_csv(data)


def wheat_leaf_db():
    from openalea.adel.wheat.extract_wheat import extract_leaf_info
    import openalea.adel.fitting as fitting

    leaves = extract_leaf_info(
        datadir + "/data/So99.RData", datadir + "/data/SRSo.RData"
    )
    leaves, discard = fitting.fit_leaves(leaves, 9)
    return leaves


def axeTable():
    filename = datadir + "/data/axeTable.csv"
    return pandas.read_csv(filename, sep=",", decimal=".")


def phytoT():
    filename = datadir + "/data/phytoT.csv"
    return pandas.read_csv(filename, sep=",", decimal=".")


def canopy_two_metamers():
    return {
        "plant": [1, 1],
        "axe_id": ["MS", "T1"],
        "refplant_id": [1, 1],
        "nff": [10, 8],
        "HS_final": [10, 8],
        "ms_insertion": [0, 1],
        "az_insertion": [0, 0],
        "numphy": [1, 1],
        "Laz": [0, 90],
        "Ll": [3, 3],
        "Lv": [3, 3],
        "Lr": [0, 0],
        "Lsen": [0, 0],
        "L_shape": [3, 3],
        "Lw_shape": [0.3, 0.3],
        "Linc": [0, 0],
        "Einc": [0, 45],
        "El": [1, 1],
        "Ev": [1, 1],
        "Esen": [0, 0],
        "Ed": [0.1, 0.1],
        "Gv": [1, 1],
        "Gl": [1, 1],
        "Gd": [0.1, 0.1],
        "Gsen": [0, 0],
        "LcType": [1, 1],
        "LcIndex": [1, 1],
    }


def adel_two_metamers(leaf_sectors=1):
    """create a very simple adel mtg"""
    ls = leaves()
    d = canopy_two_metamers()
    g = mtg_factory(d, adel_metamer, leaves=ls, leaf_sectors=leaf_sectors)
    g = mtg_interpreter(g, leaves=ls)
    return g


def adel_two_metamers_stand(
    leaf_sectors=1,
    inter_row=0.2,
    density=150,
    convunit=100,
    interleaf=1,
    leaf_length=3,
    leaf_width=0.3,
    Einc=45,
):
    """create a very simple adel mtg"""

    d = canopy_two_metamers()
    d.update(
        {
            "Ll": [leaf_length, leaf_length],
            "Lv": [leaf_length, leaf_length],
            "L_shape": [leaf_length, leaf_length],
            "Lw_shape": [leaf_width, leaf_width],
            "Einc": [0, Einc],
            "El": [1, interleaf],
            "Ev": [1, interleaf],
        }
    )

    inter_plant = 1.0 / inter_row / density
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = regular(1, 1, dx, dy)
    xc = float(domain[1][0] + domain[0][0]) / 2
    yc = float(domain[1][1] + domain[0][1]) / 2
    positions = [(x - xc, y - yc, z) for x, y, z in positions]
    domain = (
        (domain[0][0] - xc, domain[0][1] - yc),
        (domain[1][0] - xc, domain[1][1] - yc),
    )
    domain_area = (
        abs(domain[1][0] - domain[0][0])
        / convunit
        * abs(domain[1][1] - domain[0][1])
        / convunit
    )
    ls = leaves()
    g = mtg_factory(
        d, adel_metamer, leaves=ls, leaf_sectors=leaf_sectors, stand=[(positions[0], 0)]
    )
    g = mtg_interpreter(g, leaves=ls)

    return g, domain_area, domain, 1.0 / convunit


def adel_one_leaf(L=30, w=0.3, leaf_sectors=1):
    """create a very simple adel mtg"""
    ls = leaves()
    d = {
        "plant": [1],
        "axe_id": ["MS"],
        "nff": [10],
        "HS_final": [10],
        "ms_insertion": [0],
        "numphy": [1],
        "Laz": [0],
        "Ll": [3],
        "Lv": [3],
        "Lr": [0],
        "Lsen": [0],
        "L_shape": [L],
        "Lw_shape": [w],
        "Linc": [0],
        "Einc": [0],
        "El": [1],
        "Ev": [1],
        "Esen": [0],
        "Ed": [0.1],
        "Gv": [1],
        "Gl": [1],
        "Gsen": [0],
        "Gd": [0.1],
        "LcType": [1],
        "LcIndex": [1],
    }
    g = mtg_factory(d, adel_metamer, leaves=ls, leaf_sectors=leaf_sectors)
    g = mtg_interpreter(g, leaves=ls)
    return g


def adel_one_leaf_element():
    """create a very simple adel mtg"""
    ls = leaves()
    d = {
        "plant": [1],
        "axe_id": ["MS"],
        "nff": [10],
        "HS_final": [10],
        "ms_insertion": [0],
        "numphy": [1],
        "Laz": [0],
        "Ll": [3],
        "Lv": [3],
        "Lr": [0],
        "Lsen": [0],
        "L_shape": [3],
        "Lw_shape": [0.3],
        "Linc": [0],
        "Einc": [0],
        "El": [1],
        "Ev": [1],
        "Esen": [0],
        "Ed": [0.1],
        "Gv": [1],
        "Gl": [1],
        "Gsen": [0],
        "Gd": [0.1],
        "LcType": [1],
        "LcIndex": [1],
    }
    g = mtg_factory(d, adel_metamer, leaves=ls, leaf_sectors=1)
    g = mtg_interpreter(g, leaves=ls)
    g.remove_vertex(13)
    labels = g.property("label")
    labels[13] = "Removed"
    return g
