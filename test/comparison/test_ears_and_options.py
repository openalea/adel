"""ShootArchitecture's ears (peduncle, ear, awns: AD-11) and its geometry options (the blades' twist, the classic
stems, the face-up blades: AD-44) against the legacy's geometry. The legacy's static factory (newmtg.mtg_factory,
which reads Adel.R's ear rows, m_type peduncle, ear, awn, as one organ each) builds the canopy of the port's main stem
from the port's organs' lengths (in cm) and mtg_interpreter draws it; each element's triangles are compared (as
sorted vertex sets: PlantGL's tesselation orders them its own way).

Known difference (corrected): for a blade element that is its whole mature blade (fully visible, one flat part from 0
to 1), Adel's leaf_element returns the abscissae in its length unit instead of relative ones, so its twist is the
twist times the blade's length in cm; the port twists every element by the twist times its relative abscissa."""
import numpy as np
import pandas as pd
import pytest

from openalea.adel.comparison.legacy_run import triangles as plantgl_triangles
from openalea.adel.legacy import data_samples
from openalea.adel.legacy.mtg_interpreter import mtg_interpreter
from openalea.adel.legacy.newmtg import mtg_factory
from openalea.adel.simulation import simulate
from openalea.adel.tables import sample_tables, wheat_tables

CM = 100.
DEFAULTS = dict(twist=0., stem_mesh="slim", face_up=False)
OPTIONS = {"default": {}, "twist": dict(twist=30.), "classic": dict(stem_mesh="classic"),
           "face_up": dict(face_up=True), "all": dict(twist=30., stem_mesh="classic", face_up=True)}


def main_stem_run(days, inclined=False):
    """Adel's sample main stem (its ear rows: nf + 1 to nf + 3) grown *days* days at 15 °C; *inclined*: its first
    internode inclined 30° (static Adel's rules), its blades off the vertical."""
    axes, phytomers = sample_tables()
    axes = axes.loc[axes["axe"] == "MS"].reset_index(drop=True)
    phytomers = phytomers.loc[phytomers["axe"] == "MS"].reset_index(drop=True)
    if inclined:
        axes.loc[0, "incT"] = 30.
    return simulate({"axes_table": axes, "phytomer_table": phytomers, "tiller_inclination": int(inclined)},
                    days=days)


@pytest.fixture(scope="module")
def headed():
    return main_stem_run(100)


@pytest.fixture(scope="module")
def headed_inclined():
    return main_stem_run(100, inclined=True)


def canopy_table(result):
    """The legacy's canopy table (Adel.R's getdesc columns) of the port's main stem, lengths in cm."""
    ds, names = result.data_structure, result.architecture.organ_names()

    def get(name):
        return np.asarray(ds.get(name))
    by_rank = {}
    for organ, (name, rank) in enumerate(zip(names, get("organ_rank"))):
        by_rank.setdefault(int(rank), {})[name] = organ
    rows = []
    for rank in sorted(by_rank):
        organs = by_rank[rank]
        row = dict(plant=1, axe_id="MS", numphy=rank, ms_insertion=0, refplant_id=1, nff=12, HS_final=12., hasEar=1,
                   az_insertion=0., Ll=0., Lv=0., Lr=0., Lsen=0., L_shape=0., Lw_shape=0., LcType=-999, LcIndex=-999,
                   Linc=0., Laz=0., Gl=0., Gv=0., Gsen=0., Gd=0., Ginc=0., El=0., Ev=0., Esen=0., Ed=0., Einc=0.,
                   m_type="vegetative")
        stem = organs.get("internode", next((organs[k] for k in ("peduncle", "ear", "awn") if k in organs), None))
        row.update(El=get("organ_length")[stem] * CM, Ev=get("visible_length")[stem] * CM,
                   Esen=get("organ_senesced_length")[stem] * CM, Ed=get("diameter")[stem] * CM,
                   Laz=get("organ_azimuth")[stem], Einc=get("inclination")[stem])
        if "internode" in organs:
            sheath, blade = organs["sheath"], organs["blade"]
            row.update(Gl=get("organ_length")[sheath] * CM, Gv=get("visible_length")[sheath] * CM,
                       Gsen=get("organ_senesced_length")[sheath] * CM, Gd=get("diameter")[sheath] * CM,
                       Ginc=get("inclination")[sheath], Ll=get("organ_length")[blade] * CM,
                       Lv=get("visible_length")[blade] * CM, Lr=get("rolled_length")[blade] * CM,
                       Lsen=get("organ_senesced_length")[blade] * CM,
                       L_shape=get("shape_mature_length")[blade] * CM, Lw_shape=get("shape_max_width")[blade] * CM,
                       LcType=int(get("shape_class")[blade]) if get("shape_class")[blade] > 0 else -999,
                       LcIndex=int(get("shape_seed")[blade]) if get("shape_seed")[blade] > 0 else -999)
        else:
            row["m_type"] = names[stem]
        rows.append(row)
    table = pd.DataFrame(rows)
    return {column: table[column].tolist() for column in table.columns}


def port_elements(result):
    """{(rank, organ, element label): node} of the port's main stem; the plant's base (its first element's frame)."""
    ds, names = result.data_structure, result.architecture.organ_names()
    organ_of, ranks = np.asarray(ds.owner("Organ")), np.asarray(ds.get("organ_rank"))
    labels, index = np.asarray(ds.labels()), np.asarray(ds.get("element_index"))
    keys = {}
    for node, label in enumerate(labels):
        keys[(int(ranks[organ_of[node]]), names[organ_of[node]],
              label + (str(int(index[node])) if label == "LeafElement" else ""))] = node
    first = keys[(1, "internode", "HiddenElement")]
    base = np.array([np.asarray(ds.get(name))[first] for name in ("frame_x", "frame_y", "frame_z")])
    return keys, base


def compare(result, options):
    """{(rank, organ, element): max |port - legacy| (m)} of every element the legacy draws, the port's options and
    the legacy's (Leaves.twist, classic, face_up) set alike; the port's options set back after."""
    architecture = result.architecture
    leaves = data_samples.leaves()
    leaves = leaves[0] if isinstance(leaves, dict) else leaves
    leaves.twist = options.get("twist", 0.)
    g = mtg_factory(canopy_table(result), leaves={0: leaves}, leaf_sectors=1)
    g = mtg_interpreter(g, {0: leaves}, classic=options.get("stem_mesh") == "classic",
                        face_up=options.get("face_up", False))
    for name, default in DEFAULTS.items():
        setattr(architecture, name, options.get(name, default))
    try:
        triangles = architecture.triangles()
    finally:
        for name, default in DEFAULTS.items():
            setattr(architecture, name, default)
    keys, base = port_elements(result)
    differences = {}
    for vid, shape in g.property("geometry").items():
        organ = g.complex(vid)
        key = (int(g.label(g.complex(organ))[len("metamer"):]), g.label(organ), g.label(vid))
        legacy = plantgl_triangles(shape) / CM
        mine = triangles.vertices[triangles.owner == keys[key]] - base
        assert mine.shape == legacy.shape, key
        differences[key] = float(np.abs(np.sort(mine.reshape(-1, 3), axis=0)
                                        - np.sort(legacy.reshape(-1, 3), axis=0)).max()) if legacy.size else 0.
    return differences


def whole_blade(result, key):
    """Whether the element *key* is a blade element that is its whole mature blade (the legacy's twist quirk)."""
    ds = result.data_structure
    keys, _ = port_elements(result)
    node = keys[key]
    organ = np.asarray(ds.owner("Organ"))[node]
    return (key[1] == "blade" and np.isclose(np.asarray(ds.get("visible_length"))[organ],
                                             np.asarray(ds.get("shape_mature_length"))[organ], rtol=0, atol=1e-12)
            and np.asarray(ds.get("flat_base"))[node] == 0 and np.asarray(ds.get("flat_top"))[node] == 1)


# ── The ears ──

def test_the_ear_phytomers_follow_the_flag_leaf(headed):
    ds, names = headed.data_structure, headed.architecture.organ_names()
    organs = headed.organs()
    ears = organs.loc[organs["organ"].isin(["peduncle", "ear", "awn"])]
    assert ears["organ"].tolist() == ["peduncle", "ear", "awn"]
    assert ears["rank"].tolist() == [13, 14, 15]
    np.testing.assert_allclose(ears["length"], [0.198, 0.0738, 0.])            # the phytomer table's El
    ntop, diameter = np.asarray(ds.get("ntop")), np.asarray(ds.get("diameter"))
    assert ntop[names == "ear"].tolist() == [-1] and ntop[names == "awn"].tolist() == [-2]
    np.testing.assert_allclose(diameter[names == "ear"], [0.0102])
    g = ds.structure
    labels = g.labels
    peduncle = np.asarray(ds.entity_ids("Organ"))[names == "peduncle"][0]
    assert g.label(int(peduncle)) == labels.Organ.Peduncle
    hidden = g.components(int(peduncle))[0]
    flag_internode = np.asarray(ds.entity_ids("Organ"))[(names == "internode")
                                                        & (np.asarray(ds.get("organ_rank")) == 12)][0]
    assert g.parent(hidden) == g.components(int(flag_internode))[-1]        # '<' on the flag internode's top
    assert g.edge_type(hidden) == "<"


def test_the_ears_are_drawn_and_have_their_area(headed):
    ds, names = headed.data_structure, headed.architecture.organ_names()
    triangles = headed.architecture.triangles()
    organ = np.asarray(ds.owner("Organ"))
    ear_nodes = np.flatnonzero((names[organ] == "ear") & (np.asarray(ds.labels()) == "StemElement"))
    assert np.isin(ear_nodes, triangles.owner).all()
    area = np.asarray(ds.get("organ_area"))[names == "ear"][0]
    assert area == pytest.approx(np.pi * 0.0738 * 0.0102)                     # its cylinder: its whole length out


def test_no_ears_without_their_rows():
    axes, phytomers = wheat_tables(leaf_number=6, tiller_number=0)            # hasEar 1, no ear rows
    result = simulate({"axes_table": axes, "phytomer_table": phytomers}, days=60)
    assert set(result.organs()["organ"]) == {"internode", "sheath", "blade"}


def test_no_ears_on_an_axis_without_one():
    axes, phytomers = sample_tables()
    axes = axes.loc[axes["axe"] == "MS"].assign(hasEar=0).reset_index(drop=True)
    result = simulate({"axes_table": axes, "phytomer_table": phytomers.loc[phytomers["axe"] == "MS"]}, days=100)
    assert set(result.organs()["organ"]) == {"internode", "sheath", "blade"}


# ── Against the legacy ──

@pytest.mark.parametrize("option", ["default", "classic"])
def test_the_ears_and_stems_are_the_legacy_s(headed, option):
    differences = compare(headed, OPTIONS[option])
    assert {key[1] for key in differences} >= {"peduncle", "ear"}
    assert max(differences.values()) < 1e-10


@pytest.mark.parametrize("option", ["face_up", "all"])
def test_the_inclined_stem_s_options_are_the_legacy_s(headed_inclined, option):
    differences = compare(headed_inclined, OPTIONS[option])
    for key, difference in differences.items():
        if not (OPTIONS[option].get("twist") and whole_blade(headed_inclined, key)):
            assert difference < 1e-10, key


def test_the_twist_is_the_legacy_s_but_on_whole_blades():
    young = main_stem_run(60)                          # blades still growing: partly visible
    differences = compare(young, OPTIONS["twist"])
    quirk = [key for key in differences if whole_blade(young, key)]
    assert all(differences[key] < 1e-10 for key in differences if key not in quirk)
    assert any(differences[key] > 1e-6 for key in differences if key[1] == "blade")      # twisted indeed


# ── The defaults and the options ──

def test_the_options_change_what_they_name(headed_inclined):
    architecture = headed_inclined.architecture
    ds = headed_inclined.data_structure
    default = architecture.triangles()
    leaf = np.asarray(ds.labels())[default.owner] == "LeafElement"
    for name, value in (("twist", 30.), ("face_up", True)):
        setattr(architecture, name, value)
        try:
            changed = architecture.triangles()
        finally:
            setattr(architecture, name, DEFAULTS[name])
        np.testing.assert_array_equal(changed.vertices[~leaf], default.vertices[~leaf])
        assert not np.array_equal(changed.vertices[leaf], default.vertices[leaf])
    architecture.stem_mesh = "classic"
    try:
        classic = architecture.triangles()
    finally:
        architecture.stem_mesh = "slim"
    stems = np.unique(default.owner[~leaf])
    assert (np.bincount(classic.owner, minlength=ds.n_nodes())[stems] == 24).all()
    np.testing.assert_allclose(classic.areas()[np.isin(classic.owner, stems)].sum(),
                               default.areas()[np.isin(default.owner, stems)].sum(), rtol=0.3)


def test_the_options_are_set_by_the_scenario():
    result = simulate({"twist": 30., "stem_mesh": "classic", "face_up": True}, days=20)
    assert (result.architecture.twist, result.architecture.stem_mesh, result.architecture.face_up) == (30., "classic",
                                                                                                        True)
