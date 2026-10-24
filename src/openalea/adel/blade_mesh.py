"""
The blade elements' meshes, as Adel draws them (geometric_elements.Leaves.mesh: incline_leaf, then fitting.mesh4's
leaf_element and leaf_to_mesh), in the element's local frame (x along its up, z along its heading), in m.

    triangles = blade_element_triangles(shape, inclination, mature_length, max_width, visible_length, s_base, s_top,
                                        min_area)

A blade's visible part is drawn as the shape's top: its midrib from the shape's base for the visible length, its
width from the shape's width profile at the same distances from the tip (the blade grows out of the tube). An element
draws its flat part [s_base, s_top] (relative to the visible length): its midrib's points in that part, with those
nearer than (s_top - s_base) / (2 n) to the previous one dropped, the width cut at its first zero. Two rows of points
at ± half the width across the midrib make its triangles; those of area below min_area are left out.
"""
from math import atan2, cos, pi, sin

import numpy as np

TIP_WIDTH = 1e-5                            # m, Adel's 0.001 cm: a narrower tip ends the mesh with one triangle


def incline(shape, inclination):
    """Adel's incline_leaf (a relative angle): the midrib turned about the origin so that its base angle to the
    vertical becomes *inclination* times its own (at most π), then moved by its first point."""
    x, y, s, r = shape
    base_angle = atan2(x[1] - x[0], y[1] - y[0])
    rotation = base_angle - min(pi, inclination * base_angle)
    c, si = cos(rotation), sin(rotation)
    return x[0] + c * x - si * y, y[0] + si * x + c * y, s, r


def blade_part(shape, mature_length, visible_length, s_base, s_top, max_width):
    """Adel's fitting.leaf_element: the midrib's points (x, y), their abscissa and the width r of the part [s_base,
    s_top] of the visible length, scaled to the blade's mature length and maximal width; None when nothing is
    visible."""
    s_base = min(s_base, s_top, 1.)
    s_top = max(s_base, s_top, 0.)
    if visible_length <= 0.:
        return None
    visible_length = min(visible_length, mature_length)
    x, y, s, r = (np.array(values, dtype=float) for values in shape)
    if visible_length == mature_length and s_base == 0 and s_top == 1:
        return x * mature_length, y * mature_length, s * mature_length, r * max_width
    r[-1] = 0.
    visible = visible_length / mature_length
    # the midrib's abscissa from the base, the width's the same distances from the tip
    points = np.unique(np.concatenate([s, [1 - visible, visible]]))
    midrib = points[points <= visible]
    width = points[points >= 1 - visible]
    midrib = np.union1d(midrib, width - width[0])
    tip_offset = width[0]
    s_base, s_top = s_base * visible, s_top * visible
    inside = np.unique(np.concatenate([midrib, [s_base, s_top]]))
    inside = inside[(inside <= s_top) & (inside >= s_base)]
    # the points too near the previous one dropped
    tolerance = (s_top - s_base) / (len(s) * 2)
    kept = inside[1:][np.diff(inside) >= tolerance]
    kept = np.unique(np.concatenate([kept, [inside[0], inside[-1]]]))
    x_part, y_part = np.interp(kept, s, x), np.interp(kept, s, y)
    r_part = np.interp(kept + tip_offset, s, r)
    zero = r_part <= 0
    if zero.any():
        end = np.searchsorted(zero, True)               # as Adel's (a search on the mask)
        x_part, y_part, r_part = x_part[:end + 1], y_part[:end + 1], r_part[:end + 1]
        kept = kept[:end + 1]
        r_part[-1] = 0.
    return x_part * mature_length, y_part * mature_length, kept, r_part * max_width


def strip_triangles(x, y, r):
    """Adel's leaf_to_mesh without twist: points (x, ∓r/2, y) in two rows, two triangles between consecutive points,
    one at the tip when it is narrower than TIP_WIDTH; (2n, 3) points, (T, 3) indices."""
    n = len(x)
    points = np.concatenate([np.stack([x, -r / 2., y + 0.], axis=1), np.stack([x, r / 2., y - 0.], axis=1)])
    index = np.arange(n - 2) if n > 2 else np.array([0])
    triangles = [np.stack([index, index + n, index + n + 1], axis=1), np.stack([index, index + n + 1, index + 1],
                                                                              axis=1)]
    if r[-1] < TIP_WIDTH:
        triangles.append(np.array([[n - 2, 2 * n - 2, 2 * n - 1]]))
    else:
        triangles.append(np.array([[n - 2, 2 * n - 2, 2 * n - 1], [n - 2, 2 * n - 1, n - 1]]))
    return points, np.concatenate(triangles)


def blade_element_triangles(shape, inclination, mature_length, max_width, visible_length, s_base, s_top, min_area):
    """A blade element's triangles (T, 3, 3) in its local frame (Adel's Leaves.mesh with flipx): *inclination* a
    factor of the shape's base angle (a dynamic database's shapes carry their inclination: 1)."""
    x, y, s, r = incline(shape, inclination)
    part = blade_part((-x, y, s, r), mature_length, visible_length, s_base, s_top, max_width)
    if part is None or len(part[0]) < 2:
        return np.zeros((0, 3, 3))
    x_part, y_part, _, r_part = part
    points, indices = strip_triangles(x_part, y_part, r_part)
    triangles = points[indices]
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    area = np.linalg.norm(np.cross(b - a, c - a), axis=1) / 2.
    return triangles[area > min_area]
