from openalea.plantgl.all import *
from openalea.adel.symbol import *
from openalea.adel.mtg import *
import openalea.adel.fitting as fitting
import openalea.adel.json_numpy as json_np

symbols = {
    "newPlant": 1,
    "newAxe": 2,
    "newMetamer": 3,
    "StemElement": 4,
    "LeafElement": 4,
}


def leaves_db():
    import openalea.adel.fitting as fitting
    from openalea.adel import data
    from os.path import join

    pth = data.__path__[0]
    fn = join(pth, "simpleleavesdb.json")
    with open(fn) as f:
        leaves = json_np.load(f)
    leaves, discard = fitting.fit_leaves(leaves, 9)
    assert leaves is not None
    return leaves


functions = build_symbols(leaves_db())


def save_table(table, fn):
    from rpy import r

    r["write.table"](table, fn)


def save_file(str, fn):
    f = open(fn, "w")
    f.write(str)
    f.close()


def test1():
    s = """
newPlant
[newAxe
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.540000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.095000,0.308000,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.559999,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.075000,0.401579,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,1.361539,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.000000,0.522500,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.838462,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,11.800000,0.766154,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.982353,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.000000,1.015882,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.417647,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,18.900000,1.062400,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,1.300000,0.04,0.04)StemElement(1,12.999998,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.647058,1.135625,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,4.700003,0.04,0.04)StemElement(1,14.000000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,20.799999,1.235000,0.000000,1,0,0.5)]
 ]
StemElement(1,0.436111,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,7.665278,0.346528,0.000000,1,10,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(-45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.398890,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.525000,0.312500,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.588077,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.850000,0.424400,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.448352,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.750000,0.640714,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.323809,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,14.100000,0.900952,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.104762,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.500000,1.037917,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,10.650002,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,19.900000,1.116250,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,13.549996,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,22.100000,1.195714,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,3.200003,0.04,0.04)StemElement(1,16.500000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,16.500000,1.325714,0.000000,1,0,0.5)]
 ]
StemElement(1,0.000000,0.04,0.04)[/(0.000000)+(0.993730)LeafElement(1,10.050000,0.383524,0.006270,1,9,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.000000,0.04,0.04)[/(180.000000)+(0.999426)LeafElement(1,7.900000,0.455263,0.000574,1,8,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.743889,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.800000,0.358000,0.000000,1,7,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.679999,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,8.550000,0.523333,0.000000,1,6,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,2.518572,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,9.700000,0.714286,0.000000,1,5,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.344506,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,13.100000,1.014286,0.000000,1,4,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.076922,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,17.400000,1.104286,0.000000,1,3,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,10.900000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.200001,1.196552,0.000000,1,2,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,13.800000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,23.216667,1.302667,0.000000,1,1,0.5)]
newMetamer
StemElement(1,3.200001,0.04,0.04)StemElement(1,16.250000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,19.549999,1.438667,0.000000,1,0,0.5)]
]
"""
    g = CanMTG(functions, s)
    table = g.to_aggregation_table()
    # save_table(table, 'rtable.txt')
    canstr = g.to_canestra()
    # save_file(canstr, 'wheat.can')
    scene = g.to_plantgl()
    # Viewer.display(scene)
    # raw_input('enter')


def test2():
    s = """
newPlant
[newAxe
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.540000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.095000,0.308000,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.559999,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.075000,0.401579,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,1.361539,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.000000,0.522500,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.838462,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,11.800000,0.766154,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.982353,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.000000,1.015882,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.417647,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,18.900000,1.062400,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,1.300000,0.04,0.04)StemElement(1,12.999998,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.647058,1.135625,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,4.700003,0.04,0.04)StemElement(1,14.000000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,20.799999,1.235000,0.000000,1,0,0.5)]
 ]
StemElement(1,0.436111,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,7.665278,0.346528,0.000000,1,10,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(-45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.398890,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.525000,0.312500,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.588077,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.850000,0.424400,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.448352,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.750000,0.640714,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.323809,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,14.100000,0.900952,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.104762,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.500000,1.037917,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,10.650002,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,19.900000,1.116250,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,13.549996,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,22.100000,1.195714,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,3.200003,0.04,0.04)StemElement(1,16.500000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,16.500000,1.325714,0.000000,1,0,0.5)]
 ]
StemElement(1,0.000000,0.04,0.04)[/(0.000000)+(0.993730)LeafElement(1,10.050000,0.383524,0.006270,1,9,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.000000,0.04,0.04)[/(180.000000)+(0.999426)LeafElement(1,7.900000,0.455263,0.000574,1,8,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.743889,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.800000,0.358000,0.000000,1,7,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.679999,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,8.550000,0.523333,0.000000,1,6,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,2.518572,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,9.700000,0.714286,0.000000,1,5,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.344506,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,13.100000,1.014286,0.000000,1,4,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.076922,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,17.400000,1.104286,0.000000,1,3,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,10.900000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.200001,1.196552,0.000000,1,2,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,13.800000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,23.216667,1.302667,0.000000,1,1,0.5)]
newMetamer
StemElement(1,3.200001,0.04,0.04)StemElement(1,16.250000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,19.549999,1.438667,0.000000,1,0,0.5)]
]
"""
    s = s + s + s
    g = CanMTG(functions, s)
    distribution = [(0, 0, 0), (0, 90, 0), (90, 0, 0)]
    g.planter(distribution)
    scene = g.to_plantgl()
    table = g.to_aggregation_table()
    # save_table(table, 'rtable.txt')
    canstr = g.to_canestra()
    # save_file(canstr, 'wheat.can')


def test3():
    s = """
newPlant
[newAxe
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.540000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.095000,0.308000,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.559999,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.075000,0.401579,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,1.361539,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.000000,0.522500,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.838462,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,11.800000,0.766154,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.982353,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.000000,1.015882,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,6.417647,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,18.900000,1.062400,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,1.300000,0.04,0.04)StemElement(1,12.999998,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.647058,1.135625,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,4.700003,0.04,0.04)StemElement(1,14.000000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,20.799999,1.235000,0.000000,1,0,0.5)]
 ]
StemElement(1,0.436111,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,7.665278,0.346528,0.000000,1,10,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)
 [+(-45)newAxe
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.398890,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,5.525000,0.312500,0.000000,1,7,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,0.588077,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.850000,0.424400,0.000000,1,6,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.448352,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,9.750000,0.640714,0.000000,1,5,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.323809,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,14.100000,0.900952,0.000000,1,4,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,3.104762,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,17.500000,1.037917,0.000000,1,3,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,10.650002,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,19.900000,1.116250,0.000000,1,2,0.5)]
 newMetamer
 StemElement(1,0.000000,0.04,0.04)StemElement(1,13.549996,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,22.100000,1.195714,0.000000,1,1,0.5)]
 newMetamer
 StemElement(1,3.200003,0.04,0.04)StemElement(1,16.500000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,16.500000,1.325714,0.000000,1,0,0.5)]
 ]
StemElement(1,0.000000,0.04,0.04)[/(0.000000)+(0.993730)LeafElement(1,10.050000,0.383524,0.006270,1,9,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.000000,0.04,0.04)[/(180.000000)+(0.999426)LeafElement(1,7.900000,0.455263,0.000574,1,8,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.743889,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,7.800000,0.358000,0.000000,1,7,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,0.679999,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,8.550000,0.523333,0.000000,1,6,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,2.518572,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,9.700000,0.714286,0.000000,1,5,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.344506,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,13.100000,1.014286,0.000000,1,4,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,3.076922,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,17.400000,1.104286,0.000000,1,3,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,10.900000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,20.200001,1.196552,0.000000,1,2,0.5)]
newMetamer
StemElement(1,0.000000,0.04,0.04)StemElement(1,13.800000,0.04,0.04)[/(0.000000)+(1.000000)LeafElement(1,23.216667,1.302667,0.000000,1,1,0.5)]
newMetamer
StemElement(1,3.200001,0.04,0.04)StemElement(1,16.250000,0.04,0.04)[/(180.000000)+(1.000000)LeafElement(1,19.549999,1.438667,0.000000,1,0,0.5)]
]
"""
    s = s * 24
    g = CanMTG(functions, s)
    distribution = [(20 * i, 20 * j, 0) for i in range(4) for j in range(6)]
    g.planter(distribution)
    scene = g.to_plantgl()
    table = g.to_aggregation_table()
    # save_table(table, 'rtable.txt')
    canstr = g.to_canestra()
    # save_file(canstr, 'wheat.can')

    # Viewer.display(scene)
    # raw_input('enter')
