"""
The blade elements' meshes, as Adel draws them (geometric_elements.Leaves.mesh: incline_leaf, then fitting.mesh4's
leaf_element and leaf_to_mesh), in the element's local frame (x along its up, z along its heading), in m.

    triangles = blade_element_triangles(x, y, s, r, inclination, mature_length, max_width, visible_length, s_base,
                                        s_top, min_area)            the shape's fitted midrib (x, y, s, r)

A blade's visible part is drawn as the shape's top: its midrib from the shape's base for the visible length, its
width from the shape's width profile at the same distances from the tip (the blade grows out of the tube). An element
draws its flat part [s_base, s_top] (relative to the visible length): its midrib's points in that part, with those
nearer than (s_top - s_base) / (2 n) to the previous one dropped, the width cut at its first zero. Two rows of points
at ± half the width across the midrib make its triangles; those of area below min_area are left out.
"""
from math import atan2, cos, pi, sin

import numpy as np
from numba import njit

TIP_WIDTH = 1e-5                            # m, Adel's 0.001 cm: a narrower tip ends the mesh with one triangle


@njit(cache=True)
def incline(x, y, inclination):
    """Adel's incline_leaf (a relative angle): the midrib turned about the origin so that its base angle to the
    vertical becomes *inclination* times its own (at most π), then moved by its first point."""
    base_angle = atan2(x[1] - x[0], y[1] - y[0])
    rotation = base_angle - min(pi, inclination * base_angle)
    c, si = cos(rotation), sin(rotation)
    return x[0] + c * x - si * y, y[0] + si * x + c * y


@njit(cache=True)
def blade_part(x, y, s, r, mature_length, visible_length, s_base, s_top, max_width):
    """Adel's fitting.leaf_element: the midrib's points (x, y) and the width r of the part [s_base, s_top] of the
    visible length, scaled to the blade's mature length and maximal width (empty when nothing is visible)."""
    s_base = min(s_base, s_top, 1.)
    s_top = max(s_base, s_top, 0.)
    if visible_length <= 0.:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    visible_length = min(visible_length, mature_length)
    if visible_length == mature_length and s_base == 0 and s_top == 1:
        return x * mature_length, y * mature_length, r * max_width
    r = r.copy()
    r[-1] = 0.
    visible = visible_length / mature_length
    # the midrib's abscissa from the base, the width's the same distances from the tip
    points = np.unique(np.concatenate((s, np.array([1 - visible, visible]))))
    midrib = points[points <= visible]
    width = points[points >= 1 - visible]
    tip_offset = width[0]
    midrib = np.unique(np.concatenate((midrib, width - tip_offset)))
    s_base, s_top = s_base * visible, s_top * visible
    inside = np.unique(np.concatenate((midrib, np.array([s_base, s_top]))))
    inside = inside[(inside <= s_top) & (inside >= s_base)]
    # the points too near the previous one dropped
    tolerance = (s_top - s_base) / (len(s) * 2)
    kept = inside[1:][np.diff(inside) >= tolerance]
    kept = np.unique(np.concatenate((kept, np.array([inside[0], inside[-1]]))))
    x_part, y_part = np.interp(kept, s, x), np.interp(kept, s, y)
    r_part = np.interp(kept + tip_offset, s, r)
    zero = (r_part <= 0).astype(np.int64)
    if zero.any():
        end = np.searchsorted(zero, 1)                  # as Adel's (a search on the mask)
        x_part, y_part, r_part = x_part[:end + 1], y_part[:end + 1], r_part[:end + 1].copy()
        r_part[-1] = 0.
    return x_part * mature_length, y_part * mature_length, r_part * max_width


@njit(cache=True)
def strip_triangles(x, y, r):
    """Adel's leaf_to_mesh without twist: points (x, ∓r/2, y) in two rows, two triangles between consecutive points,
    one at the tip when it is narrower than TIP_WIDTH; (T, 3, 3)."""
    n = len(x)
    points = np.zeros((2 * n, 3))
    points[:n, 0], points[:n, 1], points[:n, 2] = x, -r / 2., y + 0.
    points[n:, 0], points[n:, 1], points[n:, 2] = x, r / 2., y - 0.
    pairs = max(n - 2, 1)
    tip = 1 if r[-1] < TIP_WIDTH else 2
    indices = np.zeros((2 * pairs + tip, 3), dtype=np.int64)
    for i in range(pairs):
        indices[i] = (i, i + n, i + n + 1)
        indices[pairs + i] = (i, i + n + 1, i + 1)
    indices[2 * pairs] = (n - 2, 2 * n - 2, 2 * n - 1)
    if tip == 2:
        indices[2 * pairs + 1] = (n - 2, 2 * n - 1, n - 1)
    triangles = np.zeros((indices.shape[0], 3, 3))
    for t in range(indices.shape[0]):
        for k in range(3):
            triangles[t, k] = points[indices[t, k]]
    return triangles


@njit(cache=True)
def blade_element_triangles(x, y, s, r, inclination, mature_length, max_width, visible_length, s_base, s_top,
                            min_area):
    """A blade element's triangles (T, 3, 3) in its local frame (Adel's Leaves.mesh with flipx), from its shape's
    fitted midrib (x, y, s, r): *inclination* a factor of the shape's base angle (a dynamic database's shapes carry
    their inclination: 1); those of area at most *min_area* left out."""
    x_inclined, y_inclined = incline(x, y, inclination)
    x_part, y_part, r_part = blade_part(-x_inclined, y_inclined, s, r, mature_length, visible_length, s_base, s_top,
                                        max_width)
    if len(x_part) < 2:
        return np.zeros((0, 3, 3))
    triangles = strip_triangles(x_part, y_part, r_part)
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    u, v = b - a, c - a
    area = np.sqrt((u[:, 1] * v[:, 2] - u[:, 2] * v[:, 1]) ** 2 + (u[:, 2] * v[:, 0] - u[:, 0] * v[:, 2]) ** 2
                   + (u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]) ** 2) / 2.
    return triangles[area > min_area]
