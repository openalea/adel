from openalea.plantgl import all as pgl

DISTANCE = "BOX"
DISTANCE = "MESH"
DISTANCE = "LOVEL"


def my_distance(n, sectors):
    """Compute the distance with infectable sectors."""
    if "distance" not in n._g._properties:
        n._g.add_property("distance")
    tesselator = pgl.Tesselator()
    bbc = pgl.BBoxComputer(tesselator)

    def base(geom):
        if geom:
            pts = geom.pointList
            if pts:
                n = len(pts)
                assert n % 2 == 0
                return 1 / 2 * (pts[0] + pts[n / 2])

    def top(geom):
        if geom:
            pts = geom.pointList
            if pts:
                n = len(pts)
                assert n % 2 == 0
                return 1 / 2 * (pts[n / 2 - 1] + pts[-1])

    def bbox_distance(g1, g2):
        bbc.process(pgl.Scene([g1]))
        bbox1 = bbc.result
        bbc.process(pgl.Scene([g2]))
        bbox2 = bbc.result

        return bbox1.distance(bbox2)

    def mesh_distance(g1, g2):
        pts1 = g1.pointList
        pts2 = g2.pointList
        if pts1 and pts2:
            d = maxint
            for pt in pts2:
                p1, i = pts1.findClosest(pt)
                # p1 = pts1[i]
                d = min(pgl.norm(pt - p1), d)
                return d
        else:
            return maxint

    def lovel_distance(g1, g2):
        p1 = base(g1)
        p2 = top(g2)
        return pgl.norm(p2 - p1)

    def opt_distance(n1, n2):
        g1 = n1.geometry
        g2 = n2.geometry
        if not g1 or not g2:
            return maxint
        if DISTANCE == "BOX":
            return bbox_distance(g1, g2)
        elif DISTANCE == "LOVEL":
            return lovel_distance(g1, g2)
        elif DISTANCE == "MESH":
            return mesh_distance(g1, g2)

    dists = [opt_distance(n, source) for source in sectors]
    dist = min(dists) if dists else maxint

    if not n.distance:
        n.distance = []
    if dist != maxint:
        n.distance.append((n.age, dist))

    print(n, dist)
    if dist <= 1:
        n.color = 0, 0, 255


def my_distance(n, sectors):
    """Compute the distance with infectable sectors."""
    if "distance" not in n._g._properties:
        n._g.add_property("distance")
    tesselator = pgl.Tesselator()
    bbc = pgl.BBoxComputer(tesselator)

    def base(geom):
        if geom:
            pts = geom.pointList
            if pts:
                n = len(pts)
                assert n % 2 == 0
                return 1 / 2 * (pts[0] + pts[n / 2])

    def top(geom):
        if geom:
            pts = geom.pointList
            if pts:
                n = len(pts)
                assert n % 2 == 0
                return 1 / 2 * (pts[n / 2 - 1] + pts[-1])

    def bbox_distance(g1, g2):
        bbc.process(pgl.Scene([g1]))
        bbox1 = bbc.result
        bbc.process(pgl.Scene([g2]))
        bbox2 = bbc.result

        # TODO bbox.lowerLeftCorner, bbox.upperRightCorner
        return bbox1.distance(bbox2)

    def mesh_distance(g1, g2):
        pts1 = g1.pointList
        pts2 = g2.pointList
        if pts1 and pts2:
            d = maxint
            for pt in pts2:
                p1, i = pts1.findClosest(pt)
                # TODO: To Fix
                d = pgl.norm(pt - p1)
                return d
        else:
            return maxint

    def lovel_distance(g1, g2):
        p1 = base(g1)
        p2 = top(g2)
        return pgl.norm(p2 - p1)

    def opt_distance(n1, n2):
        g1 = n1.geometry
        g2 = n2.geometry
        if not g1 or not g2:
            return maxint
        if DISTANCE == "BOX":
            return bbox_distance(g1, g2)
        elif DISTANCE == "LOVEL":
            return lovel_distance(g1, g2)
        elif DISTANCE == "MESH":
            return mesh_distance(g1, g2)

    dists = [opt_distance(n, source) for source in sectors]
    dist = min(dists) if dists else maxint

    if not n.distance:
        n.distance = []
    if dist != maxint:
        n.distance.append((n.age, dist))

    print(n, dist)
    if dist <= 1:
        n.color = 0, 0, 255
