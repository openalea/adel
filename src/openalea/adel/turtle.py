"""
The placement of ShootArchitecture's elements: Adel's turtle interpretation (mtg_interpreter.AdelVisitor), compiled,
on arrays.

The elements are visited plant by plant, each plant from its base, in Adel's order: depth first, a node's branches
before its successor (openalea.mtg's pre_order2_with_filter); traversal_order gives it, once per topology. The turtle
does not follow the topology itself: each element starts from the frame on top of the elements placed so far (the
context's top), which the visible parts move up:

- a plant starts at its position, heading up, its up towards -x, rolled by the plant's rotation;
- a tiller's first element (its axis' first phytomer, at the axis' first visit: Adel's metamer0) starts on the main
  stem's top when the previous element was the main stem's (the main stem's top and the tillers' base saved there);
  otherwise on the tillers' base; any other change of axis continues on the main stem's top saved: the main stem's
  (T1 is visited after the coleoptile: traversal_order), and a tiller's after its own tiller, a tiller of a tiller
  (AD-37); with own_axis_top (the correction), an axis continuing starts on its own top where the turtle left it (a
  tiller with its first inclined element still to come if it was), and a tiller on the top its topological parent
  left: a tiller of a tiller on its bearer's phytomer;
- a stem element with an inclination is turned down by it: an axis' first one along its azimuth (an axis' after
  another tiller rolled first by its own azimuth), the others towards the vertical plane first, then rolled back
  into their azimuth plane; a stem element with an azimuth is then rolled by it;
- after a stem element, the top is its end; after a blade element, the end of its rolled part when it has one;
  hidden elements leave the top as it is.

interpret returns each element's frame where it is drawn (its position, heading and up) and its end (the position
plus its length along the heading for a stem element, its rolled part's for a blade element). interpret and
stem_inclinations are numba kernels, in openalea.adel.kernels (the elements' kinds HIDDEN, STEM, LEAF there).

The stem organs' inclinations, by static Adel's rules (Adel.R's getIncT and axe_inclination; AdelDyn sets them to 0):
an axis' first internode inclined by its base inclination (a tiller's ramps with its Haun stage, base_inclination),
and a tiller straightened (stem_inclinations): where its stem, inclined, reaches its straightening distance from its
insertion horizontally, the stem element there is bent so that its end falls at that distance, and the next one
back to the vertical.
"""
from math import acos, cos, degrees, pi, radians

import numpy as np

from openalea.metafspm.kernels import down, left_of, roll_right, roll_to_vertical


def traversal_order(ds):
    """(order, frame_source): the positions of the DataStructure's nodes (the elements, in its order) in Adel's
    traversal: per plant, from its base element (the one without topological parent among the nodes: none, or another
    model's vertex it is anchored on), depth first, branches before successors (each group in the MPG's children
    order): ds.preorder.

    One change to Adel's: a tiller borne on the main stem's phytomer 1 (T1, in the axil of leaf 1) is visited after
    that phytomer's sheath and blade, as the last branch of its internode's last element, instead of right after the
    internode's first element it is branched on (ds.preorder's defer). The phytomer's sheath is the coleoptile, whose
    base is below ground: T1 then starts at its top, where phytomer 2 (and T2) starts, as it appears, oriented as at
    its branching point: frame_source gives its first element the position of the element it is branched on, whose
    heading and up it takes (-1 for the others). The MPG keeps T1 on phytomer 1."""
    g = ds.structure
    nodes = ds.entity_ids("node").tolist()
    position = {vid: i for i, vid in enumerate(nodes)}
    edge_type = g.topology_arrays()["edge_type"]
    branch = g._EDGE_TYPE_CODES["+"]
    roots = [vid for vid in nodes if g.parent(vid) not in position]
    frame_source = np.full(len(nodes), -1, dtype=np.int64)
    defer = []                                         # (a coleoptile tiller, the element it is visited after)
    for root in roots:
        last = g.components(g.complex(root))[-1]
        for tiller in g.children(root):
            if tiller in position and edge_type[tiller] == branch:
                defer.append((position[tiller], position[last]))
                frame_source[position[tiller]] = position[root]
    return ds.preorder(branches_first=True, defer=defer), frame_source


def base_inclination(table_inclination, is_main, haun_stage, start_haun_stage=1., rate=30.):
    """Each axis' first internode's inclination (°, Adel.R's getIncT): the table's (incT) for a main stem; a tiller's
    3° until its Haun stage reaches *start_haun_stage*, then *rate* degrees per Haun stage, from 3° to the table's."""
    ramp = np.maximum(3., np.minimum(table_inclination, rate * (haun_stage - start_haun_stage)))
    return np.where(is_main, table_inclination, np.where(haun_stage > start_haun_stage, ramp, 3.))


