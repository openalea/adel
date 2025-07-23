from heapq import *
from openalea.plantgl.all import Vector3

points = [
    Vector3(*pt)
    for pt in zip(list(range(10)), list(range(5)) + list(range(5, 0, -1)), [0] * 10)
]


def max_distance(pts, line):
    d_line = line.__normSquared__()
    max_dist = 0.0
    index = 0
    for i, pt in enumerate(pts):
        d = (Vector3(pt) ^ line).__normSquared__()
        if d > max_dist:
            max_dist = d
            index = i

    return index, max_dist / d_line


def distance(pt, p0, p1):
    line = p1 - p0
    length = line.__normSquared__()
    d = ((pt - p0) ^ line).__normSquared__()
    d /= length
    return d


def cost(polyline, nb_points):
    nb_points += 2
    n = len(polyline)
    pts = [pt for pt in polyline]
    _cost = []

    sibling = [[i - 1, i + 1] for i in range(n)]

    # compute the cost for each points
    heap_cost = []
    _cost = {}

    for i in range(1, n - 1):
        d = distance(polyline[i], polyline[i - 1], polyline[i + 1])
        _cost[i] = d
        heappush(heap_cost, [d, i])

    while len(heap_cost) > nb_points - 2:
        d, i = heappop(heap_cost)

        assert pts[i] is not None
        pts[i] = None

        # update i-1 and i+1 distance
        il, ir = sibling[i]

        if il != 0:
            sibling[il][1] = ir
            ill = sibling[il][0]
            dl = distance(pts[il], pts[ill], pts[ir])
            old_dl = _cost[il]
            heap_index = heap_cost.index([old_dl, il])
            _cost[il] = dl
            del heap_cost[heap_index]
            heappush(heap_cost, [dl, il])
        if ir != n - 1:
            sibling[ir][0] = il
            irr = sibling[ir][1]
            dr = distance(pts[ir], pts[il], pts[irr])
            old_dr = _cost[ir]
            heap_index = heap_cost.index([old_dr, ir])
            _cost[ir] = dr
            del heap_cost[heap_index]
            heappush(heap_cost, [dr, ir])
    # bug in the algorithm...
    pts[1] = None
    pts[-2] = None
    return pts
