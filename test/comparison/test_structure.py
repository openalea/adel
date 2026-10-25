"""ShootArchitecture's structure against the legacy's (AdelDyn), on the harness' scenario: every axis' (the main stem
and its tillers) phytomers and organs at each step, their elements (the legacy's base and top sentinels aside; the
legacy keeps the elements its organs' current lengths need, a subsequence of the port's stable ones), dimensions (the
legacy in cm, the port in m) and shapes."""
import numpy as np
import pytest

from openalea.adel.comparison.legacy_run import legacy_run
from openalea.adel.comparison.port_run import port_run

STEPS, EVERY = 120, 20
SENTINELS = ("baseElement", "topElement")


def chain(record, organ):
    """The organ's elements in topological order: from its first (the one whose parent is outside the organ), along
    the children inside it."""
    inside = [e for e in record.vertices if record.scale[e] == 5 and record.complex[e] == organ]
    children = {e: [c for c in inside if record.parent[c] == e] for e in inside}
    first = [e for e in inside if record.parent[e] not in inside]
    ordered = []
    while first:
        e = first.pop(0)
        ordered.append(e)
        first = children[e] + first
    return ordered


def legacy_structure(record, leaves):
    """{(rank, organ): {...}} of the legacy record's main stem, in the port's units."""
    out = {}
    for organ in record.vertices:
        if record.scale[organ] != 4 or record.label[organ] not in ("internode", "sheath", "blade"):
            continue
        metamer = record.complex[organ]
        rank = int(record.label[metamer][len("metamer"):])
        axis = record.label[record.complex[metamer]]
        elements = [record.label[e].rstrip("0123456789") for e in chain(record, organ)       # LeafElement1: its index
                    if record.label[e] not in SENTINELS]
        props = {name: values.get(organ) for name, values in record.props.items()}
        out[(axis, rank, record.label[organ])] = dict(labels=elements, diameter=(props.get("diameter") or 0.) * 0.01,
                                                shape_mature_length=(props.get("shape_mature_length") or 0.) * 0.01,
                                                shape_max_width=(props.get("shape_max_width") or 0.) * 0.01,
                                                shape_key=props.get("shape_key"), ntop=props.get("ntop"),
                                                azimuth=props.get("azimuth"))
    return out


@pytest.fixture(scope="module", params=["samples", "echap"])
def runs(request):
    leaves = request.param
    return (leaves, legacy_run(n_steps=STEPS, leaves=leaves, record_every=EVERY),
            port_run(n_steps=STEPS, leaves=leaves, record_every=EVERY))


def test_the_same_organs_and_elements(runs):
    _, legacy, port = runs
    for legacy_record, port_record in zip(legacy.steps, port.steps):
        expected = legacy_structure(legacy_record, None)
        assert set(port_record) == set(expected)
        for key, organ in expected.items():
            # the legacy rebuilds the elements of the organs' current lengths: a subsequence of the port's, stable
            port_labels = iter(port_record[key]["labels"])
            assert all(label in port_labels for label in organ["labels"]), (key, organ["labels"],
                                                                            port_record[key]["labels"])


def test_the_same_dimensions(runs):
    leaves, legacy, port = runs
    adel = legacy.adel
    keys = sorted(adel.leaves[0].srdb.keys()) if hasattr(adel.leaves[0], "srdb") else None
    expected = legacy_structure(legacy.steps[-1], leaves)
    got = port.steps[-1]
    for (axis, rank, name), organ in expected.items():
        port_organ = got[(axis, rank, name)]
        assert port_organ["ntop"] == organ["ntop"]
        if name in ("internode", "sheath"):
            assert port_organ["diameter"] == pytest.approx(organ["diameter"], rel=1e-12)
        if name == "internode":
            assert port_organ["organ_azimuth"] == pytest.approx(organ["azimuth"], rel=1e-12)
        if name == "blade":
            assert port_organ["shape_mature_length"] == pytest.approx(organ["shape_mature_length"], rel=1e-12)
            assert port_organ["shape_max_width"] == pytest.approx(organ["shape_max_width"], rel=1e-12)
            key = organ["shape_key"]
            assert (keys[port_organ["shape_class"] - 1], port_organ["shape_seed"] - 1) == tuple(key[:2])


def insertion_of_legacy(record, axis_label):
    """(bearer axis, bearer rank) of a legacy axis: up the topological parents of its first metamer's internode
    elements until an element of another axis (through the axis' metamer0 / collar)."""
    def axis_of(element):
        return record.label[record.complex[record.complex[record.complex[element]]]]

    def rank_of(element):
        return int(record.label[record.complex[record.complex[element]]][len("metamer"):])

    first = [e for e in record.vertices if record.scale[e] == 5 and axis_of(e) == axis_label and rank_of(e) == 1
             and record.label[record.complex[e]] == "internode"]
    element = min(first)
    while record.parent[element] is not None and axis_of(record.parent[element]) == axis_label:
        element = record.parent[element]
    bearer = record.parent[element]
    return axis_of(bearer), rank_of(bearer)


def test_the_tillers_on_the_same_phytomers(runs):
    _, legacy, port = runs
    record = legacy.steps[-1]
    tillers = sorted({record.label[v] for v in record.vertices if record.scale[v] == 2 and record.label[v] != "MS"})
    assert tillers
    scene = port.scene
    ds = scene.populations[0].data_structure
    g = ds._mtg
    from openalea.adel.comparison.port_run import axis_names
    names = dict(zip(np.asarray(ds.entity_ids("Axis")).tolist(), axis_names(scene)))
    for tiller in tillers:
        axis = [a for a, name in names.items() if name == tiller][0]
        first_phytomer = g.components(g.components(axis)[0])[0]
        first_element = g.components(g.components(first_phytomer)[0])[0]
        bearer_element = g.parent(first_element)
        bearer_phytomer = g.complex_at_scale(bearer_element, g.scales.Phytomer)
        bearer_axis = g.complex_at_scale(bearer_element, g.scales.Axis)
        rank = g.components(g.components(bearer_axis)[0]).index(bearer_phytomer) + 1
        assert g.edge_type(first_element) == "+"
        assert (names[bearer_axis], rank) == insertion_of_legacy(record, tiller), tiller
