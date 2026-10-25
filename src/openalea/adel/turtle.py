"""
The placement of ShootArchitecture's elements: Adel's turtle interpretation (mtg_interpreter.AdelVisitor), compiled,
on arrays.

The elements are visited plant by plant, each plant from its base, in Adel's order: depth first, a node's branches
before its successor (openalea.mtg's pre_order2_with_filter); traversal_order gives it, once per topology. The turtle
does not follow the topology itself: each element starts from the frame on top of the elements placed so far (the
context's top), which the visible parts move up:

- a plant starts at its position, heading up, its up towards -x, rolled by the plant's rotation;
- a tiller's first element (its axis' first phytomer) starts on the main stem's top when the previous element was
  the main stem's (the main stem's top and the tillers' base saved there); otherwise on the tillers' base; the
  main stem continues on the main stem's top saved (T1 is visited after the coleoptile: traversal_order);
- a stem element with an inclination is turned down by it: an axis' first one along its azimuth (an axis' after
  another tiller rolled first by its own azimuth), the others towards the vertical plane first, then rolled back
  into their azimuth plane; a stem element with an azimuth is then rolled by it;
- after a stem element, the top is its end; after a blade element, the end of its rolled part when it has one;
  hidden elements leave the top as it is.

interpret returns each element's frame where it is drawn (its position, heading and up) and its end (the position
plus its length along the heading for a stem element, its rolled part's for a blade element).

The stem organs' inclinations, by static Adel's rules (Adel.R's getIncT and axe_inclination; AdelDyn sets them to 0):
an axis' first internode inclined by its base inclination (a tiller's ramps with its Haun stage, base_inclination),
and a tiller straightened (stem_inclinations): where its stem, inclined, reaches its straightening distance from its
insertion horizontally, the stem element there is bent so that its end falls at that distance, and the next one
back to the vertical.
"""
from math import acos, cos, degrees, pi, radians

import numpy as np
from numba import njit

from openalea.metafspm.coupling.geometry import down, left_of, roll_right, roll_to_vertical

HIDDEN, STEM, LEAF = 0, 1, 2


def traversal_order(g, nodes):
    """(order, frame_source): the positions of *nodes* (the elements' vertex ids, in the DataStructure's order) in
    Adel's traversal: per
    plant, from its base element (the one without topological parent), depth first, branches before successors (each
    group in the MPG's children order).

    One change to Adel's: a tiller borne on the main stem's phytomer 1 (T1, in the axil of leaf 1) is visited after
    that phytomer's sheath and blade, as the last branch of its internode's last element, instead of right after the
    internode's first element it is branched on. The phytomer's sheath is the coleoptile, whose base is below ground:
    T1 then starts at its top, where phytomer 2 (and T2) starts, as it appears, oriented as at its branching point:
    frame_source gives its first element the position of the element it is branched on, whose heading and up it
    takes (-1 for the others). The MPG keeps T1 on phytomer 1."""
    position = {vid: i for i, vid in enumerate(nodes)}
    edge_type = g.property("edge_type")
    roots = [vid for vid in nodes if g.parent(vid) is None]
    deferred = {}                                      # an element: the coleoptile tillers visited after its branches
    for root in roots:
        tillers = [child for child in g.children(root) if child in position and edge_type.get(child) == "+"]
        if tillers:
            deferred[g.components(g.complex(root))[-1]] = tillers
    moved = {tiller for tillers in deferred.values() for tiller in tillers}
    frame_source = np.full(len(nodes), -1, dtype=np.int64)
    for root in roots:
        for tiller in deferred.get(g.components(g.complex(root))[-1], []):
            frame_source[position[tiller]] = position[root]
    order = []
    for root in roots:
        stack = [root]
        while stack:
            vid = stack.pop()
            order.append(position[vid])
            children = [child for child in g.children(vid) if child in position and child not in moved]
            branches = [child for child in children if edge_type.get(child) != "<"] + deferred.get(vid, [])
            successors = [child for child in children if edge_type.get(child) == "<"]
            stack.extend(reversed(branches + successors))
    return np.asarray(order, dtype=np.int64), frame_source


@njit(cache=True)
def interpret(order, kind, plant_base, axis, axis_is_main, tiller_start, plant_position, plant_rotation,
              axis_azimuth, inclination, azimuth, length, rolled, frame_source):
    """Per element: (position, heading, up, end), each (n, 3), by Adel's turtle (above). *order*: traversal_order;
    per element: *kind* (HIDDEN, STEM, LEAF), *plant_base* (its plant's first), *axis* (its axis' index),
    *tiller_start* (on a tiller's first phytomer), its plant's *plant_position* (n, 3) and *plant_rotation* (°), its
    axis' *axis_azimuth* (°), its *inclination* and *azimuth* (°, stem elements), *length*, *rolled* (a blade
    element's rolled part), *frame_source* (traversal_order's: a deferred tiller's heading and up are that element's);
    per axis: *axis_is_main*."""
    n = kind.size
    position, heading, up, end = np.zeros((n, 3)), np.zeros((n, 3)), np.zeros((n, 3)), np.zeros((n, 3))
    top = np.zeros((3, 3))                       # rows: position, heading, up
    main_top, tiller_base = np.zeros((3, 3)), np.zeros((3, 3))
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
            is_axis_first = True
        if a != previous_axis:
            if tiller_start[i]:
                is_axis_first = True
                if axis_is_main[previous_axis]:
                    main_top[:] = top
                    if frame_source[i] >= 0:           # at the coleoptile's top, oriented as where it branches
                        top[1], top[2] = heading[frame_source[i]], up[frame_source[i]]
                    tiller_base[:] = top
                else:
                    top[:] = tiller_base
            else:
                top[:] = main_top
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
        previous_axis = a
    return position, heading, up, end


def base_inclination(table_inclination, is_main, haun_stage, start_haun_stage=1., rate=30.):
    """Each axis' first internode's inclination (°, Adel.R's getIncT): the table's (incT) for a main stem; a tiller's
    3° until its Haun stage reaches *start_haun_stage*, then *rate* degrees per Haun stage, from 3° to the table's."""
    ramp = np.maximum(3., np.minimum(table_inclination, rate * (haun_stage - start_haun_stage)))
    return np.where(is_main, table_inclination, np.where(haun_stage > start_haun_stage, ramp, 3.))


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
