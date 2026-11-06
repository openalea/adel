"""
Adel's numba kernels, in one place (metafspm's conventions.md, "Performance rules"): the turtle's interpretation of
the stems and the blade meshes, each @njit(cache=True). ShootArchitecture calls them (openalea.adel.turtle
describes the turtle). warm_up() runs every one on tiny inputs, which fills numba's on-disk cache (metafspm's
warm_up runs it with every port's: python -m openalea.metafspm warm_up).
"""
from math import acos, atan2, cos, degrees, pi, radians, sin

import numpy as np
from numba import njit

from openalea.metafspm.kernels import down, left_of, roll_right, roll_to_vertical

HIDDEN, STEM, LEAF = 0, 1, 2                # the elements' kinds (interpret)
TIP_WIDTH = 1e-5                            # m, Adel's 0.001 cm: a narrower tip ends the mesh with one triangle


# ── The turtle (openalea.adel.turtle): the elements' frames, the stems' inclinations ──


@njit(cache=True)
def interpret(order, kind, plant_base, axis, axis_is_main, tiller_start, plant_position, plant_rotation,
              axis_azimuth, inclination, azimuth, length, rolled, frame_source, parent, own_axis_top):
    """Per element: (position, heading, up, end), each (n, 3), by Adel's turtle (above). *order*: traversal_order;
    per element: *kind* (HIDDEN, STEM, LEAF), *plant_base* (its plant's first), *axis* (its axis' index),
    *tiller_start* (on a tiller's first phytomer), its plant's *plant_position* (n, 3) and *plant_rotation* (°), its
    axis' *axis_azimuth* (°), its *inclination* and *azimuth* (°, stem elements), *length*, *rolled* (a blade
    element's rolled part), *frame_source* (traversal_order's: a deferred tiller's heading and up are that element's),
    *parent* (its topological parent's position, -1 for none), *own_axis_top* (its plant's correction AD-37: at a
    change of axis, an axis continuing on its own top, a tiller on the top its parent left); per axis:
    *axis_is_main*."""
    n = kind.size
    position, heading, up, end = np.zeros((n, 3)), np.zeros((n, 3)), np.zeros((n, 3)), np.zeros((n, 3))
    left_top = np.zeros((n, 3, 3))               # per element: the top it left
    top = np.zeros((3, 3))                       # rows: position, heading, up
    main_top, tiller_base = np.zeros((3, 3)), np.zeros((3, 3))
    tiller_continuation = np.zeros((3, 3))       # the legacy's main stem's top saved for a tiller after its tiller
    entered = np.zeros(axis_is_main.size, dtype=np.bool_)
    left_axis_top = np.zeros((axis_is_main.size, 3, 3))           # per axis: the top when the turtle left it
    left_first = np.zeros(axis_is_main.size, dtype=np.bool_)      # and is_axis_first
    is_axis_first = True
    previous_axis = -1
    for i in order:
        a = axis[i]
        if previous_axis < 0:
            previous_axis = a
        if plant_base[i]:
            top[0] = plant_position[i]
            top[1] = np.array([0., 0., 1.])
            top[2] = np.array([-1., 0., 0.])
            roll_right(top[1], top[2], plant_rotation[i])
            main_top[:] = top
            tiller_base[:] = top
            tiller_continuation[:] = top
            is_axis_first = True
        if a != previous_axis:
            left_axis_top[previous_axis], left_first[previous_axis] = top, is_axis_first
            if tiller_start[i] and not entered[a]:     # the legacy's metamer0: a tiller's first visit only
                is_axis_first = True
                if axis_is_main[previous_axis]:
                    main_top[:] = top
                    if frame_source[i] >= 0:           # at the coleoptile's top, oriented as where it branches
                        top[1], top[2] = heading[frame_source[i]], up[frame_source[i]]
                    tiller_base[:] = top
                    tiller_continuation[:] = top           # main_top but for a deferred tiller's orientation
                else:
                    top[:] = tiller_base
            elif axis_is_main[a]:
                top[:] = main_top
            else:
                top[:] = tiller_continuation
            if own_axis_top[i] and entered[a]:                  # an axis continuing: as the turtle left it
                top[:] = left_axis_top[a]
                if not axis_is_main[a]:
                    is_axis_first = left_first[a]
            elif own_axis_top[i] and parent[i] >= 0 and frame_source[i] < 0:     # a tiller: on its bearer's top
                top[:] = left_top[parent[i]]
        entered[a] = True
        p, h, u = top[0].copy(), top[1].copy(), top[2].copy()
        if kind[i] == STEM:
            if inclination[i] != 0.:
                if is_axis_first:
                    if not axis_is_main[previous_axis]:
                        roll_right(h, u, axis_azimuth[i])
                        tiller_base[0], tiller_base[1], tiller_base[2] = p, h, u
                    down(h, u, inclination[i])
                    is_axis_first = False
                else:
                    up_before = u.copy()
                    left_z = left_of(h, u)[2]
                    roll_to_vertical(h, u)
                    cosine = min(1., max(-1., up_before[0] * u[0] + up_before[1] * u[1] + up_before[2] * u[2]))
                    angle = degrees(acos(cosine))
                    left_change = left_z - left_of(h, u)[2]
                    down(h, u, inclination[i])
                    if left_change < 0:
                        angle = -angle
                    roll_right(h, u, -angle)
            if azimuth[i] != 0.:
                roll_right(h, u, azimuth[i])
        position[i], heading[i], up[i] = p, h, u
        if kind[i] == STEM:
            if length[i] > 0:
                p = p + h * length[i]
            top[0], top[1], top[2] = p, h, u
        elif kind[i] == LEAF:
            if rolled[i] > 0:
                p = p + h * rolled[i]
                top[0], top[1], top[2] = p, h, u
        end[i] = p
        left_top[i] = top
        previous_axis = a
    return position, heading, up, end


@njit(cache=True)
def stem_inclinations(axis, kind, visible, first_inclination, straightening, epsilon=1e-6):
    """The stem organs' inclinations (°, Adel.R's axe_inclination). *axis*, *kind* (0 internode, 1 sheath), *visible*
    (their visible lengths): per organ, sorted by axis, rank and kind (each axis from its first internode); per axis:
    its first internode's inclination and its straightening distance (Adel's dredT, in the lengths' unit; 0: none).
    The stem's visible organs follow each other from the insertion; the first whose end is at the straightening
    distance horizontally (the stem straight, at its first internode's inclination) is bent so that its end falls at
    it, and the next one is turned back to the vertical."""
    n = axis.size
    out = np.zeros(n)
    start = 0
    while start < n:
        a = axis[start]
        stop = start
        total = 0.
        while stop < n and axis[stop] == a:
            total += visible[stop]
            stop += 1
        out[start] = first_inclination[a]
        if straightening[a] > 0 and total > epsilon:
            alpha = radians(90. - first_inclination[a])
            elements = np.array([k for k in range(start, stop) if visible[k] > 0])
            length, distance = 0., 0.
            for j in range(elements.size):
                k = elements[j]
                previous_length, previous_distance = length, distance
                length += visible[k]
                distance = length * cos(alpha)
                if distance >= straightening[a]:
                    if j > 0:
                        beta = acos((straightening[a] - previous_distance) / (length - previous_length))
                        out[k] = -degrees(beta - alpha)
                        if j + 1 < elements.size:
                            out[elements[j + 1]] = -degrees(pi / 2 - beta)
                    else:                                  # the base inclination too large: the first one bent
                        beta = acos(straightening[a] / length)
                        out[start] = degrees(pi / 2 - beta)
                        if elements.size > 1:
                            out[elements[1]] = -out[start]
                    break
        start = stop
    return out


# ── The blade meshes ──
#
# The blade elements' meshes, as Adel draws them (geometric_elements.Leaves.mesh: incline_leaf, then fitting.mesh4's
# leaf_element and leaf_to_mesh), in the element's local frame (x along its up, z along its heading), in m.
#
#     triangles = blade_element_triangles(x, y, s, r, inclination, mature_length, max_width, visible_length, s_base,
#                                         s_top, min_area)            the shape's fitted midrib (x, y, s, r)
#     ends = blade_element_midrib_ends(x, y, s, r, inclination, mature_length, max_width, visible_length, s_base,
#                                      s_top)                         its midrib's first and last points
#
# A blade's visible part is drawn as the shape's top: its midrib from the shape's base for the visible length, its
# width from the shape's width profile at the same distances from the tip (the blade grows out of the tube). An element
# draws its flat part [s_base, s_top] (relative to the visible length): its midrib's points in that part, with those
# nearer than (s_top - s_base) / (2 n) to the previous one dropped, the width cut at its first zero. Two rows of points
# at ± half the width across the midrib make its triangles; those of area below min_area are left out.


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


@njit(cache=True)
def blade_element_midrib_ends(x, y, s, r, inclination, mature_length, max_width, visible_length, s_base, s_top):
    """A blade element's flat part's midrib ends (2, 3), its first and last points, in its local frame as
    blade_element_triangles draws it (its triangles' rows lie on either side of them); (0, 3) when it has no
    triangles."""
    x_inclined, y_inclined = incline(x, y, inclination)
    x_part, y_part, _ = blade_part(-x_inclined, y_inclined, s, r, mature_length, visible_length, s_base, s_top,
                                   max_width)
    ends = np.zeros((0, 3))
    if len(x_part) >= 2:
        ends = np.zeros((2, 3))
        ends[0, 0], ends[0, 2] = x_part[0], y_part[0]
        ends[1, 0], ends[1, 2] = x_part[-1], y_part[-1]
    return ends


# ── Warm-up ──

def warm_up() -> None:
    """Every kernel compiled (or loaded from numba's cache) on a tiny shoot: a short simulation of Adel's (its
    ShootArchitecture placing and meshing the elements: interpret, stem_inclinations, the blade meshes), which runs
    them on the types they are given. metafspm's warm_up runs it (python -m openalea.metafspm warm_up)."""
    from openalea.adel.simulation import simulate
    simulate(days=20)
