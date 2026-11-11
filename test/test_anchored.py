"""ShootArchitecture on a plant another model shares: a crown initiator (a stand-in for the use case's: a seed and its
subcrown internode, StemElements up to the crown) builds first, the main stem's first element follows the crown's top
(the plant's scenario's shoot_anchor), and a model of the seed works on it. ShootArchitecture sees its own elements
only (SHOOT_ELEMENTS, its `entities`), its development its own axes and organs, and the shoot grows as a plain plant
sown at the crown's depth; the crown and the seed keep their own values."""
from dataclasses import dataclass, replace

import numpy as np
import pytest

from openalea.adel.development import DAY, ScenarioShoot, shoot_scene
from openalea.metafspm.scene.population import stand
from openalea.metafspm.coupling.composite_wrapper import CompositeModel
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
    def initiate_plant(cls, mtg, plant, scenario):
        scales, labels = mtg.scales, mtg.labels
        parent = plant
        for scale, label in ((scales.Axis, labels.Axis.Shoot), (scales.GrowthUnit, labels.GrowthUnit.Shoot),
                             (scales.Phytomer, labels.Phytomer.Shoot), (scales.Organ, labels.Organ.StemInternode)):
            parent = mtg.add_component(parent, **PropsConfig(scale=scale, edge_type='/', label=label))
        x, y, z = (float(mtg.property(name)[plant]) for name in ("x", "y", "z"))
        element = mtg.add_component(parent, **PropsConfig(scale=scales.SubOrgan, edge_type='/',
                                                        label=labels.SubOrgan.Seed))
        cls._place(mtg, element, x, y, z, z)
        crown = -float(scenario.get("crown_depth", CROWN_DEPTH))
        rise = (crown - z) / SUBCROWN_ELEMENTS
        for k in range(SUBCROWN_ELEMENTS):
            element = mtg.add_child(element, **PropsConfig(scale=scales.SubOrgan, edge_type='<',
                                                         label=labels.SubOrgan.StemElement))
            cls._place(mtg, element, x, y, z + k * rise, z + (k + 1) * rise)
            mtg.set_values("length", [element], [rise])
            mtg.set_values("area", [element], [CROWN_VALUES["area"]])

    @staticmethod
    def _place(mtg, element, x, y, bottom, top):
        for name, value in dict(x1=x, y1=y, z1=bottom, x2=x, y2=y, z2=top).items():
            mtg.set_values(name, [element], [value])


@dataclass
class SeedReserve(FunctionalComponent):
    """A model of the seed only (its entities): its reserve used at a constant rate."""
    entities = {"label": "Seed"}

    reserve: float = state_variable(unit="g", description="The seed's reserve", default=0.04, kind="extensive",
                                    scale=ScalesConfig.SubOrgan)
    use_rate: float = parameter(unit="g s-1", description="The reserve's use", by="model_seed", default=1e-9)

    @state
    def _reserve(self, reserve, use_rate):
        return np.maximum(reserve - self.step_duration * use_rate, 0.)


class CrownedShoot(CompositeModel):
    """ScenarioShoot (nested) on a crown, its seed's reserve used."""
    initiators = (StandInCrown, ScenarioShoot)
    parts = {"shoot": ScenarioShoot, "seed": SeedReserve}


def _scene(plants, model, **parameters):
    return shoot_scene(plants, model=model, time_step=DAY, depth=None, forcings={"air_temperature": 15.},
                       start=0., delays=[0., 40.][:len(plants)], tiller_inclination=1, n_sect=2, **parameters)


def _grown(n_days=70):
    """(the crowned scene, the plain one): two plants sown at SOWING_DEPTH on a crown at CROWN_DEPTH, and the same
    plants without a crown, sown at CROWN_DEPTH; grown n_days."""
    plants = stand(2, rows=4, plant_spacing=0.04, depth=SOWING_DEPTH, jitter=0.005, seed=0)
    crowned = _scene(plants, CrownedShoot, shoot_anchor=CROWN_TOP, crown_depth=CROWN_DEPTH)
    plain = _scene([replace(plant, depth=CROWN_DEPTH) for plant in plants], ScenarioShoot)
    for _ in range(n_days):
        crowned.run()
        plain.run()
    return crowned, plain


@pytest.fixture(scope="module")
def grown():
    return _grown()


def _by_creation(data_structure, values, location="node"):
    """*values* (one per entity of *location*) in the entities' creation order (their vids'): the crowned and the plain
    DataStructures order their nodes along their topology, from the seed or from the main stem's base."""
    return np.asarray(values)[np.argsort(np.asarray(data_structure.entity_ids(location)), kind="stable")]


def _architecture(scene):
    """The scene's ShootArchitecture: its ScenarioShoot's (the plain model's own, the crowned one's nested)."""
    model = scene.populations[0].instance
    return getattr(model, "shoot", model).architecture


def _triangles_by_creation(scene):
    """The architecture's triangles' vertices, by their elements' creation order (each element's in their order)."""
    architecture = _architecture(scene)
    triangles = architecture.triangles()
    ids = np.asarray(architecture.data_structure.entity_ids("node"))
    rank = np.argsort(np.argsort(ids, kind="stable"), kind="stable")
    return triangles.vertices[np.argsort(rank[triangles.owner], kind="stable")]


def _crown(data_structure):
    """The crown's elements' indices among the DataStructure's nodes, its seed's, its axis' and its organ's."""
    labels = np.asarray(data_structure.labels())
    crown_axis = np.flatnonzero(np.asarray(data_structure.get("axis_row")) < 0)
    subcrown = np.flatnonzero((labels == "StemElement") & np.isin(data_structure.owner("Axis"), crown_axis))
    organ = np.unique(np.asarray(data_structure.owner("Organ"))[subcrown])
    return subcrown, np.flatnonzero(labels == "Seed"), crown_axis, organ


def test_its_elements_only_its_main_stem_on_the_crown(grown):
    crowned, _ = grown
    data_structure = crowned.populations[0].data_structure
    architecture = _architecture(crowned)
    view = architecture.data_structure
    assert isinstance(view, EntityView)
    subcrown, seed, crown_axis, _ = _crown(data_structure)
    assert len(subcrown) == 2 * SUBCROWN_ELEMENTS and len(seed) == 2 and len(crown_axis) == 2
    own = np.asarray(view.entity_ids("node"))
    subcrown_ids = np.asarray(data_structure.entity_ids("node"))[subcrown]
    assert not np.isin(np.asarray(data_structure.entity_ids("node"))[np.r_[subcrown, seed]], own).any()
    assert view.n_nodes() == data_structure.n_nodes() - 3 * 2
    mtg = data_structure.structure
    for plant, base in mtg.property("shoot_base").items():
        top = mtg.parent(base)                                           # the crown's top: no StemElement of it above
        assert top in subcrown_ids and not np.isin(mtg.children(top), subcrown_ids).any()
        assert mtg.edge_type(base) == '<' and base in own
        assert mtg.complex(mtg.complex(mtg.complex(mtg.complex(base)))) != mtg.complex(mtg.complex(mtg.complex(mtg.complex(top))))
    assert (view.parents() == -1).sum() == 2                           # each main stem's base: on its crown
    assert set(architecture.axis_names()[crown_axis].tolist()) == {""}
    assert {"MS", "T1", "T2"} <= set(architecture.axis_names().tolist())   # tillers grown


def test_it_grows_as_a_plain_plant_sown_at_the_crown(grown):
    crowned, plain = grown
    data_structure, plain_ds = crowned.populations[0].data_structure, plain.populations[0].data_structure
    view = _architecture(crowned).data_structure
    assert view.n_nodes() == plain_ds.n_nodes()
    for name in ELEMENT_VARIABLES:
        np.testing.assert_array_equal(_by_creation(view, view.get(name)), _by_creation(plain_ds, plain_ds.get(name)),
                                      err_msg=name)
    for name in PLACEMENT_VARIABLES:                 # the crown's top: its subcrown's z2, the same height rounded
        np.testing.assert_allclose(_by_creation(view, view.get(name)), _by_creation(plain_ds, plain_ds.get(name)),
                                   rtol=0., atol=1e-12, err_msg=name)
    own = {"Axis": np.asarray(data_structure.get("axis_row")) >= 0}
    own["Organ"] = own["Axis"][data_structure.owner("Axis", of="Organ")]
    for location, names in (("Organ", ORGAN_VARIABLES), ("Axis", AXIS_VARIABLES)):
        ids = np.asarray(data_structure.entity_ids(location))[own[location]]
        for name in names:
            mine = np.asarray(data_structure.get(name))[own[location]][np.argsort(ids, kind="stable")]
            np.testing.assert_array_equal(mine, _by_creation(plain_ds, plain_ds.get(name), location), err_msg=name)
    np.testing.assert_allclose(_triangles_by_creation(crowned), _triangles_by_creation(plain), rtol=0., atol=1e-12)
    assert np.isclose(np.min(np.asarray(view.get("z1"))), -CROWN_DEPTH)  # its base on the crown's top


def test_the_crown_and_the_seed_keep_their_own(grown):
    crowned, _ = grown
    data_structure = crowned.populations[0].data_structure
    subcrown, seed, crown_axis, crown_organ = _crown(data_structure)
    for name, value in CROWN_VALUES.items():
        np.testing.assert_allclose(np.asarray(data_structure.get(name))[subcrown], value, err_msg=name)
    mtg = data_structure.structure
    tops = data_structure.index_of(np.array([mtg.parent(base) for base in mtg.property("shoot_base").values()]))
    np.testing.assert_allclose(np.asarray(data_structure.get("z2"))[tops], -CROWN_DEPTH)
    np.testing.assert_array_equal(np.asarray(data_structure.get("phytomer_number"))[crown_axis], 1.)
    np.testing.assert_array_equal(np.asarray(data_structure.get("haun_stage"))[crown_axis], 0.)
    np.testing.assert_array_equal(np.asarray(data_structure.get("axis_phytomers"))[crown_axis], 0)
    for name in ("organ_length", "visible_length", "organ_area"):
        np.testing.assert_array_equal(np.asarray(data_structure.get(name))[crown_organ], 0., err_msg=name)
    reserve = crowned.populations[0].instance.seed.data_structure.get("reserve")
    np.testing.assert_allclose(reserve, 0.04 - 70 * DAY * 1e-9)        # as the seed alone
    assert np.asarray(data_structure.get("reserve"))[seed].tolist() == np.asarray(reserve).tolist()


def test_without_its_anchor_built_first_it_raises():
    class Reversed(CrownedShoot):
        initiators = (ScenarioShoot, StandInCrown)
    with pytest.raises(AnchorError, match="StandInCrown"):
        _scene(stand(1, plant_spacing=0.04, jitter=0.005, seed=0), Reversed, shoot_anchor=CROWN_TOP)
