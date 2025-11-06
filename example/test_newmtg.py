from openalea.adel.newmtg import *
from openalea.adel.mtg_interpreter import *
import openalea.plantgl.all as pgl
from openalea.adel.AdelR import devCsv, setAdel, RunAdel, genGeoLeaf, genGeoAxe
import openalea.adel.fitting as fitting

import openalea.adel.data as data
from pathlib import Path


def adelR(nplants, dd):
    devT = devCsv(
        "./data/axeTCa0N.csv",
        "./data/dimTCa0N.csv",
        "./data/phenTCa0N.csv",
        "./data/earTCa0N.csv",
        "./data/ssi2sen.csv",
    )
    geoLeaf = genGeoLeaf()
    geoAxe = genGeoAxe()
    pars = setAdel(devT, geoLeaf, geoAxe, nplants)
    cantable = RunAdel(dd, pars)
    return pars, cantable


def leaves_db(dlevel=9):
    import pickle as Pickle

    d = data.__path__[0]
    fn = r"%s/leaves_simple.db" % str(d)
    f = open(fn)
    leaves = Pickle.load(f)
    f.close()
    leaves, discard = fitting.fit_leaves(leaves, dlevel)
    return leaves


# from openalea.mtg import mtg2axialtree
# affichage : aller voir __str__dans mtg.py
# dp = {'plant':[1,1,1,1,1], 'axe' : [0,0,1,1,2], 'numphy': [1,2,1,2,1], 'Laz': [0,180,90,180,180],'Einc':[0,0,45,-45,45], 'El': [1,1,1,1,1],'Ev':[1,1,1,1,1], 'Esen':[0,0,0,0,0], 'Ed':[0.1,0.1,0.1,0.1,0.1]}
# dp = {'plant':[1,1,1,1,1], 'axe' : [0,0,1,1,2], 'numphy': [1,2,1,2,1], 'Laz': [0,180,90,180,0],'Einc':[0,0,45,-30,45], 'El': [1,1,1,1,1],'Ev':[1,1,1,1,1], 'Esen':[0,0,0,0,0], 'Ed':[0.1,0.1,0.1,0.1,0.1]}
# g=mtg_factory(dp,adel_metamer)
# g.display()

# more test
p, d = adelR(1, 1000)
leaves = leaves_db()
g = mtg_factory(
    d,
    adel_metamer,
    leaves=leaves,
    stand=[((0, 0, 0), 0), ((10, 0, 0), 90), ((0, 10, 0), 0)],
)
g = mtg_interpreter(g)
scene = plot3d(g)
pgl.Viewer.display(scene)


#
def test_leaf_shape():
    """test leaf positioning/ inclination"""
    d = {
        "plant": [1, 1],
        "axe_id": ["MS", "T1"],
        "ms_insertion": [0, 1],
        "numphy": [1, 1],
        "Laz": [0, 90],
        "Ll": [3, 3],
        "Lv": [3, 3],
        "Lsen": [0, 0],
        "L_shape": [3, 3],
        "Lw_shape": [0.3, 0.3],
        "Linc": [0, 0],
        "Einc": [0, 45],
        "El": [1, 1],
        "Ev": [1, 1],
        "Esen": [0, 0],
        "Ed": [0.1, 0.1],
    }
    g = mtg_factory(d, adel_metamer, leaf_db=leaves_db())
    g = mtg_interpreter(g)
    scene = plot3d(g)
    pgl.Viewer.display(scene)
    return g
