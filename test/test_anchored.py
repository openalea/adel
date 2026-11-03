"""ShootArchitecture on a plant another model shares: a crown initiator (a stand-in for the use case's: a seed and its
subcrown internode, StemElements up to the crown) builds first, the main stem's first element follows the crown's top
(the plant's scenario's shoot_anchor), and a model of the seed works on it. ShootArchitecture sees its own elements
only (SHOOT_ELEMENTS, its `entities`), its development its own axes and organs, and the shoot grows as a plain plant
sown at the crown's depth; the crown and the seed keep their own values."""
from dataclasses import dataclass, replace

import numpy as np
import pytest

from openalea.adel.architecture import ShootArchitecture
from openalea.adel.development import DAY, ScenarioShoot, shoot_scene, sown_stand
from openalea.metafspm.coupling.component import FunctionalComponent, StructuralComponent, parameter, state_variable
from openalea.metafspm.data_structure.configs import PropsConfig, ScalesConfig
from openalea.metafspm.data_structure.data_api import EntityView
from openalea.metafspm.data_structure.mpg import AnchorError
from openalea.metafspm.solve.decorator import state

SOWING_DEPTH, CROWN_DEPTH, SUBCROWN_ELEMENTS = 0.04, 0.02, 2
CROWN_TOP = {"label": "StemElement", "last": True}                  # the crown's top: the last of its StemElements
ELEMENT_VARIABLES = ("length", "area", "green_length", "green_area", "senesced_area", "rolled_part", "width",
                     "is_green")
PLACEMENT_VARIABLES = ("x1", "y1", "z1", "x2", "y2", "z2", "heading_x", "heading_y", "heading_z", "up_x", "up_y",
                       "up_z")
ORGAN_VARIABLES = ("organ_length", "visible_length", "organ_senesced_length", "organ_area", "inclination", "age")
AXIS_VARIABLES = ("phytomer_number", "haun_stage", "axis_phytomers", "first_inclination")
CROWN_VALUES = {"length": (SOWING_DEPTH - CROWN_DEPTH) / SUBCROWN_ELEMENTS, "area": 1e-5}


class StandInCrown(StructuralComponent):
    """A crown initiator: under its Axis, GrowthUnit, Phytomer and Organ (a StemInternode, the subcrown internode),
    a seed at the plant's position, then StemElements up to the crown (scenario's crown_depth), placed (x1 ... z2)
    with their lengths and areas."""

    @classmethod
    def initiate_plant(cls, g, plant, scenario):
        scales, labels = g.scales, g.labels
        parent = plant
        for scale, label in ((scales.Axis, labels.Axis.Shoot), (scales.GrowthUnit, labels.GrowthUnit.Shoot),
                             (scales.Phytomer, labels.Phytomer.Shoot), (scales.Organ, labels.Organ.StemInternode)):
            parent = g.add_component(parent, **PropsConfig(scale=scale, edge_type='/', label=label))
        x, y, z = (float(g.property(name)[plant]) for name in ("x", "y", "z"))
        element = g.add_component(parent, **PropsConfig(scale=scales.SubOrgan, edge_type='/',
                                                        label=labels.SubOrgan.Seed))
        cls._place(g, element, x, y, z, z)
        crown = -float(scenario.get("crown_depth", CROWN_DEPTH))
        rise = (crown - z) / SUBCROWN_ELEMENTS
        for k in range(SUBCROWN_ELEMENTS):
            element = g.add_child(element, **PropsConfig(scale=scales.SubOrgan, edge_type='<',
                                                         label=labels.SubOrgan.StemElement))
            cls._place(g, element, x, y, z + k * rise, z + (k + 1) * rise)
            g.property("length")[element] = rise
            g.property("area")[element] = CROWN_VALUES["area"]

    @staticmethod
    def _place(g, element, x, y, bottom, top):
        for name, value in dict(x1=x, y1=y, z1=bottom, x2=x, y2=y, z2=top).items():
            g.property(name)[element] = value


@dataclass
class SeedReserve(FunctionalComponent):
    """A model of the seed only (its entities): its reserve used at a constant rate."""
    entities = {"label": "Seed"}

    reserve: float = state_variable(unit="g", description="The seed's reserve", default=0.04, kind="extensive",
                                    scale=ScalesConfig.SubOrgan)
    use_rate: float = parameter(unit="g s-1", description="The reserve's use", by="model_seed", default=1e-9)

    @state
    def _reserve(self, reserve, use_rate):
        return np.maximum(reserve - self.dt * use_rate, 0.)


class CrownedShoot(ScenarioShoot):
    """ScenarioShoot on a crown, its seed's reserve used."""
    initiators = (StandInCrown, ShootArchitecture)

    def __init__(self, data_structure, time_step):
        super().__init__(data_structure, time_step)
        self.seed = SeedReserve(data_structure=data_structure, time_step=time_step)

    def run(self):
        super().run()
        self.seed()


def _scene(plants, model, **parameters):
    return shoot_scene(plants, model=model, time_step=DAY, sowing_depth=None, forcings={"air_temperature": 15.},
                       start=0., delays=[0., 40.][:len(plants)], tiller_inclination=1, n_sect=2, **parameters)


def _grown(n_days=70):
    """(the crowned scene, the plain one): two plants sown at SOWING_DEPTH on a crown at CROWN_DEPTH, and the same
    plants without a crown, sown at CROWN_DEPTH; grown n_days."""
    plants = [replace(plant, z=-SOWING_DEPTH) for plant in sown_stand(2)]
    crowned = _scene(plants, CrownedShoot, shoot_anchor=CROWN_TOP, crown_depth=CROWN_DEPTH)
    plain = _scene([replace(plant, z=-CROWN_DEPTH) for plant in plants], ScenarioShoot)
    for _ in range(n_days):
        crowned.run()
        plain.run()
    return crowned, plain


@pytest.fixture(scope="module")
def grown():
    return _grown()


def _by_creation(ds, values, location="node"):
    """*values* (one per entity of *location*) in the entities' creation order (their vids'): the crowned and the plain
    DataStructures order their nodes along their topology, from the seed or from the main stem's base."""
    return np.asarray(values)[np.argsort(np.asarray(ds.entity_ids(location)), kind="stable")]


def _triangles_by_creation(scene):
    """The architecture's triangles' vertices, by their elements' creation order (each element's in their order)."""
    architecture = scene.populations[0].instance.architecture
    triangles = architecture.triangles()
    ids = np.asarray(architecture.data_structure.entity_ids("node"))
    rank = np.argsort(np.argsort(ids, kind="stable"), kind="stable")
    return triangles.vertices[np.argsort(rank[triangles.owner], kind="stable")]


def _crown(ds):
    """The crown's elements' indices among the DataStructure's nodes, its seed's, its axis' and its organ's."""
    labels = np.asarray(ds.labels())
    crown_axis = np.flatnonzero(np.asarray(ds.get("axis_row")) < 0)
    subcrown = np.flatnonzero((labels == "StemElement") & np.isin(ds.owner("Axis"), crown_axis))
    organ = np.unique(np.asarray(ds.owner("Organ"))[subcrown])
    return subcrown, np.flatnonzero(labels == "Seed"), crown_axis, organ


def test_its_elements_only_its_main_stem_on_the_crown(grown):
    crowned, _ = grown
    ds = crowned.populations[0].data_structure
    architecture = crowned.populations[0].instance.architecture
    view = architecture.data_structure
    assert isinstance(view, EntityView)
    subcrown, seed, crown_axis, _ = _crown(ds)
    assert len(subcrown) == 2 * SUBCROWN_ELEMENTS and len(seed) == 2 and len(crown_axis) == 2
    own = np.asarray(view.entity_ids("node"))
    subcrown_ids = np.asarray(ds.entity_ids("node"))[subcrown]
    assert not np.isin(np.asarray(ds.entity_ids("node"))[np.r_[subcrown, seed]], own).any()
    assert view.n_nodes() == ds.n_nodes() - 3 * 2
    g = ds.structure
    for plant, base in g.property("shoot_base").items():
        top = g.parent(base)                                           # the crown's top: no StemElement of it above
        assert top in subcrown_ids and not np.isin(g.children(top), subcrown_ids).any()
        assert g.edge_type(base) == '<' and base in own
        assert g.complex(g.complex(g.complex(g.complex(base)))) != g.complex(g.complex(g.complex(g.complex(top))))
    assert (view.parents() == -1).sum() == 2                           # each main stem's base: on its crown
    assert set(architecture.axis_names()[crown_axis].tolist()) == {""}
    assert {"MS", "T1", "T2"} <= set(architecture.axis_names().tolist())   # tillers grown


def test_it_grows_as_a_plain_plant_sown_at_the_crown(grown):
    crowned, plain = grown
    ds, plain_ds = crowned.populations[0].data_structure, plain.populations[0].data_structure
    view = crowned.populations[0].instance.architecture.data_structure
    assert view.n_nodes() == plain_ds.n_nodes()
    for name in ELEMENT_VARIABLES:
        np.testing.assert_array_equal(_by_creation(view, view.get(name)), _by_creation(plain_ds, plain_ds.get(name)),
                                      err_msg=name)
    for name in PLACEMENT_VARIABLES:                 # the crown's top: its subcrown's z2, the same height rounded
        np.testing.assert_allclose(_by_creation(view, view.get(name)), _by_creation(plain_ds, plain_ds.get(name)),
                                   rtol=0., atol=1e-12, err_msg=name)
    own = {"Axis": np.asarray(ds.get("axis_row")) >= 0}
    own["Organ"] = own["Axis"][ds.owner("Axis", of="Organ")]
    for location, names in (("Organ", ORGAN_VARIABLES), ("Axis", AXIS_VARIABLES)):
        ids = np.asarray(ds.entity_ids(location))[own[location]]
        for name in names:
            mine = np.asarray(ds.get(name))[own[location]][np.argsort(ids, kind="stable")]
            np.testing.assert_array_equal(mine, _by_creation(plain_ds, plain_ds.get(name), location), err_msg=name)
    np.testing.assert_allclose(_triangles_by_creation(crowned), _triangles_by_creation(plain), rtol=0., atol=1e-12)
    assert np.isclose(np.min(np.asarray(view.get("z1"))), -CROWN_DEPTH)  # its base on the crown's top


def test_the_crown_and_the_seed_keep_their_own(grown):
    crowned, _ = grown
    ds = crowned.populations[0].data_structure
    subcrown, seed, crown_axis, crown_organ = _crown(ds)
    for name, value in CROWN_VALUES.items():
        np.testing.assert_allclose(np.asarray(ds.get(name))[subcrown], value, err_msg=name)
    g = ds.structure
    tops = ds.index_of(np.array([g.parent(base) for base in g.property("shoot_base").values()]))
    np.testing.assert_allclose(np.asarray(ds.get("z2"))[tops], -CROWN_DEPTH)
    np.testing.assert_array_equal(np.asarray(ds.get("phytomer_number"))[crown_axis], 1.)
    np.testing.assert_array_equal(np.asarray(ds.get("haun_stage"))[crown_axis], 0.)
    np.testing.assert_array_equal(np.asarray(ds.get("axis_phytomers"))[crown_axis], 0)
    for name in ("organ_length", "visible_length", "organ_area"):
        np.testing.assert_array_equal(np.asarray(ds.get(name))[crown_organ], 0., err_msg=name)
    reserve = crowned.populations[0].instance.seed.data_structure.get("reserve")
    np.testing.assert_allclose(reserve, 0.04 - 70 * DAY * 1e-9)        # as the seed alone
    assert np.asarray(ds.get("reserve"))[seed].tolist() == np.asarray(reserve).tolist()


def test_without_its_anchor_built_first_it_raises():
    class Reversed(CrownedShoot):
        initiators = (ShootArchitecture, StandInCrown)
    with pytest.raises(AnchorError, match="StandInCrown"):
        _scene([sown_stand(1)[0]], Reversed, shoot_anchor=CROWN_TOP)
