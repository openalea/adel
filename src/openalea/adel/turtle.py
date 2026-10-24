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
  main stem continues on the main stem's top saved;
- a stem element with an inclination is turned down by it: an axis' first one along its azimuth (an axis' after
  another tiller rolled first by its own azimuth), the others towards the vertical plane first, then rolled back
  into their azimuth plane; a stem element with an azimuth is then rolled by it;
- after a stem element, the top is its end; after a blade element, the end of its rolled part when it has one;
  hidden elements leave the top as it is.

interpret returns each element's frame where it is drawn (its position, heading and up) and its end (the position
plus its length along the heading for a stem element, its rolled part's for a blade element).
"""
from math import acos, degrees

import numpy as np
from numba import njit

from openalea.metafspm.coupling.geometry import down, left_of, roll_right, roll_to_vertical

HIDDEN, STEM, LEAF = 0, 1, 2


def traversal_order(g, nodes):
    """The positions of *nodes* (the elements' vertex ids, in the DataStructure's order) in Adel's traversal: per
    plant, from its base element (the one without topological parent), depth first, branches before successors (each
    group in the MPG's children order)."""
    position = {vid: i for i, vid in enumerate(nodes)}
    edge_type = g.property("edge_type")
    roots = [vid for vid in nodes if g.parent(vid) is None]
    order = []
    for root in roots:
        stack = [root]
        while stack:
            vid = stack.pop()
            order.append(position[vid])
            children = [child for child in g.children(vid) if child in position]
            branches = [child for child in children if edge_type.get(child) != "<"]
            successors = [child for child in children if edge_type.get(child) == "<"]
            stack.extend(reversed(branches + successors))
    return np.asarray(order, dtype=np.int64)


@njit(cache=True)
def interpret(order, kind, plant_base, axis, axis_is_main, tiller_start, plant_position, plant_rotation,
              axis_azimuth, inclination, azimuth, length, rolled):
    """Per element: (position, heading, up, end), each (n, 3), by Adel's turtle (above). *order*: traversal_order;
    per element: *kind* (HIDDEN, STEM, LEAF), *plant_base* (its plant's first), *axis* (its axis' index),
    *tiller_start* (on a tiller's first phytomer), its plant's *plant_position* (n, 3) and *plant_rotation* (°), its
    axis' *axis_azimuth* (°), its *inclination* and *azimuth* (°, stem elements), *length*, *rolled* (a blade
    element's rolled part); per axis: *axis_is_main*."""
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
