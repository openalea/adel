# from openalea.plantgl import all as pgl
import numpy as np
from openalea.plantgl import all as pgl


def curve_discretizer(curve, n=100):
    """ """

    u = np.linspace(curve.firstKnot, curve.lastKnot, n)
    pts = list(map(curve.getPointAt, u))
    x, y = list(zip(*pts))

    # return outputs
    return (
        x,
        y,
    )


def points2curve(x, y):
    points = list(zip(x, y))
    return pgl.NurbsCurve2D.fit(points)
