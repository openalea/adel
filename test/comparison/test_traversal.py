"""turtle.traversal_order (metafspm's data_structure.preorder with the coleoptile tillers deferred) against the per-vertex walk it
replaced, on the harness' plants with tillers (T1 included) and on two plants, at several steps of their growth: the same
positions in the same order, and the same frame sources."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import StandPlant
from openalea.adel.comparison.port_run import port_scene
from openalea.adel.turtle import traversal_order


def _walk(mtg, nodes):
    """The former traversal: per plant from its base, depth first by an explicit stack, branches before successors,
    a coleoptile tiller (on a base element, edge '+') visited after its base's organ's last element."""
    position = {vid: i for i, vid in enumerate(nodes)}
    edge_type = mtg.column("edge_type")[0]
    roots = [vid for vid in nodes if mtg.parent(vid) not in position]
    deferred = {}
    for root in roots:
        tillers = [child for child in mtg.children(root) if child in position and edge_type[child] == "+"]
        if tillers:
            deferred[mtg.components(mtg.complex(root))[-1]] = tillers
    moved = {tiller for tillers in deferred.values() for tiller in tillers}
    frame_source = np.full(len(nodes), -1, dtype=np.int64)
    for root in roots:
        for tiller in deferred.get(mtg.components(mtg.complex(root))[-1], []):
            frame_source[position[tiller]] = position[root]
    order = []
    for root in roots:
        stack = [root]
        while stack:
            vid = stack.pop()
            order.append(position[vid])
            children = [child for child in mtg.children(vid) if child in position and child not in moved]
            branches = [child for child in children if edge_type[child] != "<"] + deferred.get(vid, [])
            successors = [child for child in children if edge_type[child] == "<"]
            stack.extend(reversed(branches + successors))
    return np.asarray(order, dtype=np.int64), frame_source


@pytest.mark.parametrize("plants", [None, [StandPlant(azimuth=0.), StandPlant(x=10., azimuth=30.)]])
def test_the_traversal_is_the_walk_it_replaced(plants):
    scene = port_scene("samples", plants=plants)
    checked = 0
    for step in range(121):
        scene.run()
        if step % 15:
            continue
        data_structure = scene.populations[0].instance.architecture.data_structure
        order, frame_source = traversal_order(data_structure)
        expected_order, expected_source = _walk(data_structure.structure, data_structure.entity_ids("node").tolist())
        np.testing.assert_array_equal(order, expected_order)
        np.testing.assert_array_equal(frame_source, expected_source)
        checked += int((frame_source >= 0).any())
    assert checked > 0                                     # T1 visited after the coleoptile at some steps
