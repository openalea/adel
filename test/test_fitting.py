import openalea.adel.fitting as fitting
import openalea.adel.mtg as CanMTG
from openalea.adel.symbol import build_symbols
from numpy import compress, unique, union1d, interp
import random
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
    return leaves


db = leaves_db()
functions = build_symbols(db)


def test1():
    leaf_rank = 1
    seed = 1
    total_length = 7
    length = 7
    s_base = 0.5
    s_top = 0.6
    radius_max = 1

    rank_max = max(list(map(int, list(db.keys()))))
    rank = leaf_rank
    rank = min(rank, rank_max)
    # choisi la liste de leaves du rang, ou rang + 1 si clef absente ourag -1 ou liste vide sinon
    leaves = db.get(str(rank), db.get(str(rank + 1), db.get(str(rank - 1), [])))
    n = len(leaves)
    # if n == 0:
    #    raise "Not enough leaves data at rank %d"%rank
    random.seed(seed)
    i = random.randint(0, n - 1)
    # leaf = leaves[i]
    leaf = leaves[0]
    pts, ind = fitting.mesh4(leaf, total_length, length, s_base, s_top, radius_max)

    # x, y, r = fitting.debugmesh4(leaf, total_length, length, s_base, s_top, radius_max)

    def insert_values(a, values):
        l = a.tolist()
        l.extend(values)
        return unique(l)

    length_max = total_length
    s_base = min(s_base, s_top, 1.0)
    s_top = max(s_base, s_top, 0.0)
    x, y, s, r = leaf
    # force width to zero at the top
    r[-1] = 0
    # 1. compute s_xy and s_r for length vs length_max
    if length <= 0.0:
        sys.exit()
    if length > length_max:
        length = length_max

    param = float(length) / length_max
    n = len(s)
    # add param in s_xy ICI!!!!
    sp = insert_values(s, [1 - param, param])
    s_xy = compress(sp <= param, sp)
    s_r = compress(sp >= (1 - param), sp)
    s_xy = union1d(s_xy, s_r - s_r[0])
    s_r = s_xy + s_r[0]
    # add s_base and s_top in s

    s_base *= param
    s_top *= param
    s_valid = insert_values(s_xy, [s_base, s_top])
    s_valid = compress(s_valid <= s_top, s_valid)
    s_valid = compress(s_valid >= s_base, s_valid)

    # delete small intervals COMIT from  Here !
    eps = 1.0 / (len(s) * 2)
    ds = s_valid[1:] - s_valid[:-1]
    error = ds >= eps
    s_val = compress(error, s_valid[1:])
    s_val = insert_values(s_val, [s_valid[0], s_valid[-1]])

    # find param and 1-param in s...
    # build xf, yf, rf with the same nb of points

    s_xyf = s_val
    s_rf = s_val + s_r[0]

    xf = interp(s_xyf, s, x)
    yf = interp(s_xyf, s, y)
    rf = interp(s_rf, s, r)

    cond = rf <= 0
    if cond.any():
        # delete points with negative radius
        index = cond.searchsorted(True)
        xf = xf[: index + 1]
        yf = yf[: index + 1]
        rf = rf[: index + 1]
        rf[-1] = 0.0

        xf *= length_max
        yf *= length_max
        rf *= radius_max

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
    g = CanMTG.CanMTG(functions, s)
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
    g = CanMTG.CanMTG(functions, s)
    distribution = [(20 * i, 20 * j, 0) for i in range(4) for j in range(6)]
    g.planter(distribution)
    scene = g.to_plantgl()
    table = g.to_aggregation_table()
    canstr = g.to_canestra()

    # Viewer.display(scene)
    # raw_input('enter')
