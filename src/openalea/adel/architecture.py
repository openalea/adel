"""
Adel's shoot architecture as a metafspm StructuralComponent: the plants' topology in the MPG (axes, phytomers, organs,
elements), the organs' dimensions and shapes, their elements, and their geometry for light models and display. The
organs' lengths are given by other models (development, senescence), as cnwgrass gives them to Adel's AdelDyn.

    table = planting([Planting()], MyModel, {"axes_table": axeT, "phytomer_table": phytoT, "table_plant": 1})

with MyModel.initiators = (ShootArchitecture,). Without tables, a plant is generated from a few scalars when any is
given (its parameters leaf_number, tiller_number, longest_blade_length, longest_blade_width), else it is Adel's sample
plant (openalea.adel.tables):

    scenario={"leaf_number": 11, "tiller_number": 2, "longest_blade_length": 0.2}

The plant's structure, in the MPG:

    Plant
     └─ Axis (ShootAxis)              the main stem, each tiller (one per row of the plant's axes table)
         └─ GrowthUnit (ShootGrowthUnit)
             └─ Phytomer (ShootPhytomer)   one per rank, created when the axis' phytomer_number asks for it
                 ├─ Organ StemInternode    a HiddenElement, then a StemElement (<)
                 ├─ Organ Sheath           the same, its first element branched (+) on the internode's last
                 └─ Organ Blade            a HiddenElement (its part in the sheaths' tube), then n_sect
                                           LeafElements (<), on the sheath's last element

Axes appear when the plant's axis_number asks for them, in its axes table's order, each with its first phytomer;
phytomers when their axis' phytomer_number asks for them. A phytomer's internode follows the previous one's (<); a
tiller's first internode is branched (+) on its bearer's phytomer at its insertion rank, at its internode's base (its
first element), as Adel inserts it (Adel's names: 'T1' on the main stem's phytomer 1, 'T1.2' on T1's phytomer 2). The
elements are stable: created with their organ, their lengths and areas updated, never rebuilt.

The tables are Adel's (its axeT and phytoT, one row per axis and per axis and rank, column 'plant'; `table_plant`, a
per-plant parameter, says which of their plants a plant is). Each plant's tables are kept on the MPG at its initiation;
the components index the population's (tables.population_tables: every plant's rows at an offset; tables()).

The elements are placed by Adel's turtle (openalea.adel.turtle), from each plant's position and rotation: the frame
they are drawn in (frame_x…, heading_x…, up_x…), where they are drawn (their display ends display_x1 … display_z2: a
blade element's along its blade, its rolled part's base, its flat part's midrib's end; a hidden element's a point on
the visible parts below it), and their topological ends x1, y1, z1, x2, y2, z2 (metafspm's convention: each element's
base at its parent's end, the main stem's at the plant's base; its end its base plus its chord, a stem or hidden
element's heading times its length, a blade element's drawn chord), at each step after their lengths. triangles()
gives their triangles (metafspm's coupling.geometry.Triangles) for light models and display.
The axes are vertical, as AdelDyn's; with tiller_inclination=1, inclined by static Adel's rules from the axes table
(incT ramping with the axes' Haun stage, the straightening at dredT, the azimuths azTb and azT: _first_inclination per
axis, _inclination per organ). The turtle is Adel's, which places a tiller of a tiller (T1.1) and its bearer after it
from saved frames of other axes (AD-37); with own_axis_top=1, each axis continues on its own top and a tiller starts
on its bearer's.

The coupling contract, what other models write and read (lengths in m, on the DataStructure):

- a development model writes, per Plant, `axis_number` (its axes emerged, in its axes table's order); per Axis,
  `phytomer_number` (its phytomers initiated, 1 when it appears); per Organ, `organ_length`, `visible_length` (out of
  the sheaths' tube), a blade's `rolled_length` and `age`, and optionally a blade's `final_blade_length` and
  `final_blade_width` (its final dimensions, as elongwheat's lamina_Lmax and leaf_Wmax; NaN: not given);
- ShootArchitecture's blades' shapes are `shape_mature_length` and `shape_max_width`: the phytomer table's Ll and Lw
  at their creation, the development's final dimensions once it gives them;
- a senescence model writes, per Organ, `organ_senesced_length`;
- ShootArchitecture derives the elements from them at each call (their lengths, areas, green and senesced parts,
  placement): the elements are its own, the others read them, they do not write them;
- a growth or metabolism model reads the elements' `area`, `green_area`, `senesced_area`, the organs' `organ_area`;
- a light model reads triangles() (each triangle's element, its node index) and gives its results back per element
  (Triangles.per_owner).

On a plant other models share (a crown model's seed and subcrown internode, a root model's roots), its entities are
its elements (SHOOT_ELEMENTS: metafspm's EntityView), and the plant's scenario entry shoot_anchor (a selector for
metafspm's MPG.find_anchor) hangs the main stem's first element '<' on another initiator's vertex, its base on that
vertex's top (shoot_base_height):

    scenario={"shoot_anchor": {"label": "StemElement", "last": True}}  # initiators (Crown, ShootArchitecture): its top

Its steps per axis and organ leave another model's axes (without a row in the axes table, axis_row -1) and organs as
they are.

The organs a call creates get their lengths at the next call; a development model giving them at once (as cnwgrass
after add_metamer) gives them in its postsegmentation step, run together with ShootArchitecture's
(CompositeModel.run_together: openalea.adel.development's ScenarioShoot).
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from openalea.metafspm.coupling.component import (NewComplexes, NewNodes, StructuralComponent, display_end,
                                                   input_variable, parameter, planting_variable, setting,
                                                   state_variable)
from openalea.metafspm.coupling.geometry import HalfSpace, Triangles, place, slim_cylinders
from openalea.metafspm.data_structure.configs import PropsConfig, ScalesConfig
from openalea.metafspm.solve.decorator import postsegmentation, segmentation
from openalea.adel.kernels import blade_element_midrib_ends, blade_element_triangles, interpret, stem_inclinations
from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.tables import (phytomer_dimensions, plant_rows, plant_tables, population_tables,
                                  record_plant_tables)
from openalea.adel.turtle import base_inclination, traversal_order

BY = "model_shoot_architecture"
DESCRIPTOR = "descriptor"
SCALES = ScalesConfig
ORGANS = ("internode", "sheath", "blade")
ELEMENTS = ("HiddenElement", "StemElement", "LeafElement")


# ShootArchitecture's entities (its class attribute `entities`: metafspm's EntityView, when other models share the
# plant), in metafspm's filters' syntax: the shoot's elements (ELEMENTS) on its axes, those with a row in the axes
# table. Another model's elements are not, even of the same labels (a crown's subcrown internode, its StemElements:
# on the crown's axis, axis_row -1); on an anchored plant, they are those its main stem's first element bears.
SHOOT_ELEMENTS = {"label": ELEMENTS, "axis_row": ">=0"}


@lru_cache(maxsize=1)
def default_leaves():
    """Adel's sample leaf shapes (LeafShapes.samples()), fitted once."""
    return LeafShapes.samples()


def bearer_of(axis_name):
    """(the bearer axis' name, the insertion rank) of an axis by Adel's names: MS has none; 'T1' is on MS's phytomer
    1, 'T1.2' on T1's phytomer 2."""
    if axis_name == "MS":
        return None, 0
    code = axis_name[1:].split(".")
    return ("MS" if len(code) == 1 else "T" + ".".join(code[:-1])), int(code[-1])


def phytomer_layout(labels, n_sect):
    """A phytomer's complexes and elements, in creation order:
    complexes [(anchor, scale, label)] with anchor "unit" (the axis' GrowthUnit) or the index of an earlier complex;
    elements [(complex index, parent, edge, label, organ, index in organ)] with parent "previous" (the internode's
    predecessor: the previous internode's last element, the bearer's internode's first, or none) or the index of an
    earlier element."""
    complexes = [("unit", SCALES.Phytomer, labels.Phytomer.Shoot),
                 (0, SCALES.Organ, labels.Organ.StemInternode),
                 (0, SCALES.Organ, labels.Organ.Sheath),
                 (0, SCALES.Organ, labels.Organ.Blade)]
    elements = [(1, "previous", None, labels.SubOrgan.HiddenElement, 0, 0),
                (1, 0, "<", labels.SubOrgan.StemElement, 0, 1),
                (2, 1, "+", labels.SubOrgan.HiddenElement, 1, 0),
                (2, 2, "<", labels.SubOrgan.StemElement, 1, 1),
                (3, 3, "<", labels.SubOrgan.HiddenElement, 2, 0)]
    for k in range(n_sect):
        elements.append((3, len(elements) - 1, "<", labels.SubOrgan.LeafElement, 2, k + 1))
    return complexes, elements


@dataclass
class ShootArchitecture(StructuralComponent):
    """Adel's shoot architecture: topology, organs' dimensions and shapes, elements, geometry."""

    # Its elements, when other models share the plant (a crown's, the roots'): its view of the DataStructure
    entities = SHOOT_ELEMENTS

    # --- Inputs: what development and senescence models write ---
    axis_number: float = input_variable(unit='adim',
                                        description="The number of the plant's axes emerged, in its axes table's order "
                                                    "(the main stem first)", by='model_development', default=1.,
                                        scale=SCALES.Plant)
    phytomer_number: float = input_variable(unit='adim',
                                            description="The number of phytomers the axis has (its development)",
                                            value_comment="1 at its creation: an axis appears with its first phytomer",
                                            by='model_development', default=1., scale=SCALES.Axis)
    organ_length: float = input_variable(unit='m', description="The organ's length", by='model_development',
                                         default=0., scale=SCALES.Organ)
    visible_length: float = input_variable(unit='m', description="The organ's length out of the sheaths' tube",
                                           by='model_development', default=0., scale=SCALES.Organ)
    rolled_length: float = input_variable(unit='m', description="The blade's rolled length in the sheaths' tube",
                                          by='model_development', default=0., scale=SCALES.Organ)
    organ_senesced_length: float = input_variable(unit='m', description="The organ's senesced length",
                                                  by='model_senescence', default=0., scale=SCALES.Organ)
    age: float = input_variable(unit='°Cd', description="The blade's age (its shape class when the database has ages)",
                                by='model_development', default=0., scale=SCALES.Organ)
    haun_stage: float = input_variable(unit='adim',
                                       description="The axis' Haun stage (its tillers' inclination ramps with it)",
                                       by='model_development', default=0., scale=SCALES.Axis)
    final_blade_length: float = input_variable(unit='m', description="The blade's final length, given by development",
                                               value_comment="NaN: not given (the phytomer table's Ll)",
                                               by='model_development', default=np.nan, scale=SCALES.Organ)
    final_blade_width: float = input_variable(unit='m',
                                              description="The blade's maximal width, given by development",
                                              value_comment="NaN: not given (the phytomer table's Lw)",
                                              by='model_development', default=np.nan, scale=SCALES.Organ)
    shoot_base_height: float = input_variable(unit='m',
                                              description="The height of the main stem's base, on the top of the "
                                                          "vertex it is anchored on (its scenario's shoot_anchor)",
                                              value_comment="NaN: the plant's z (the planting's); written at "
                                                            "initiation from the anchor's z2 when its model gave it",
                                              by='model_anchor', default=np.nan, scale=SCALES.Plant, optional=True)

    # --- The axes (from the axes table) ---
    shoot_tables: object = state_variable(unit='', description="The plant's tables: (axes table, phytomer table, its "
                                          "number in them), from its scenario (tables.plant_tables)",
                                          kind=DESCRIPTOR, default=None, scale=SCALES.Plant, dtype='object')
    axis_row: int = state_variable(unit='adim', description="The axis' row in the axes table",
                                   kind=DESCRIPTOR, default=-1, scale=SCALES.Axis, dtype='int')
    bearer_row: int = state_variable(unit='adim', description="The bearer axis' row in the axes table",
                                     value_comment="-1: the main stem", kind=DESCRIPTOR, default=-1,
                                     scale=SCALES.Axis, dtype='int')
    insertion_rank: int = state_variable(unit='adim', description="The rank of the bearer's phytomer the axis is on",
                                         kind=DESCRIPTOR, default=0, scale=SCALES.Axis, dtype='int')
    axis_phytomers: int = state_variable(unit='adim', description="The number of phytomers the axis has in the MPG",
                                         kind=DESCRIPTOR, default=0, scale=SCALES.Axis, dtype='int')
    final_leaf_number: int = state_variable(unit='adim', description="The axis' final number of leaves (Adel's nf)",
                                            kind=DESCRIPTOR, default=0, scale=SCALES.Axis,
                                            dtype='int')
    final_haun_stage: float = state_variable(unit='adim', description="The axis' final Haun stage (Adel's HS_final)",
                                             kind=DESCRIPTOR, default=0., scale=SCALES.Axis)
    has_ear: int = state_variable(unit='adim', description="1 when the axis bears an ear",
                                  kind=DESCRIPTOR, default=0, scale=SCALES.Axis, dtype='int')
    axis_azimuth: float = state_variable(unit='°',
                                         description="The axis' azimuth at its insertion (Adel's az_insertion), used "
                                                     "when its first stem element is inclined",
                                         kind=DESCRIPTOR, default=0., scale=SCALES.Axis)
    table_inclination: float = state_variable(unit='°',
                                              description="The axis' first internode's inclination (Adel's incT)",
                                              value_comment="used with tiller_inclination",
                                              kind=DESCRIPTOR, default=0., scale=SCALES.Axis)
    straightening_distance: float = state_variable(unit='m',
                                                   description="The horizontal distance from its insertion at which "
                                                               "the axis straightens (Adel's dredT)",
                                                   value_comment="0: none; used with tiller_inclination",
                                                   kind=DESCRIPTOR, default=0., scale=SCALES.Axis)
    insertion_azimuth: float = state_variable(unit='°', description="The axis' insertion azimuth (Adel's azTb)",
                                              value_comment="used with tiller_inclination",
                                              kind=DESCRIPTOR, default=0., scale=SCALES.Axis)
    first_azimuth: float = state_variable(unit='°', description="The azimuth of the axis' first phytomer (Adel's azT)",
                                          value_comment="used with tiller_inclination", kind=DESCRIPTOR,
                                          default=0., scale=SCALES.Axis)
    first_inclination: float = state_variable(unit='°',
                                              description="The axis' first internode's inclination (Adel.R's "
                                                          "getIncT: incT, a tiller's ramping with its Haun stage)",
                                              value_comment="used with tiller_inclination",
                                              kind=DESCRIPTOR, default=0., scale=SCALES.Axis)

    # --- The phytomers and organs (from the phytomer table, at their creation) ---
    phytomer_rank: int = state_variable(unit='adim', description="The phytomer's rank on its axis, from 1",
                                        kind=DESCRIPTOR, default=0, scale=SCALES.Phytomer,
                                        dtype='int')
    organ_rank: int = state_variable(unit='adim', description="The organ's phytomer rank",
                                     kind=DESCRIPTOR, default=0, scale=SCALES.Organ, dtype='int')
    ntop: int = state_variable(unit='adim', description="The organ's rank from the axis' top leaf (1: the flag leaf)",
                               kind=DESCRIPTOR, default=0, scale=SCALES.Organ, dtype='int')
    diameter: float = state_variable(unit='m', description="The stem organ's diameter (phytomer table: Ed, Gd)",
                                     kind=DESCRIPTOR, default=0., scale=SCALES.Organ)
    organ_azimuth: float = state_variable(unit='°',
                                          description="The organ's azimuth (the internode's: phytomer table's Azim)",
                                          kind=DESCRIPTOR, default=0., scale=SCALES.Organ)
    inclination: float = state_variable(unit='°', description="The organ's inclination",
                                        kind=DESCRIPTOR, default=0., scale=SCALES.Organ)
    shape_class: int = state_variable(unit='adim', description="The blade's shape class in the leaf database (Lindex)",
                                      value_comment="-1: none", kind=DESCRIPTOR, default=-1,
                                      scale=SCALES.Organ, dtype='int')
    shape_seed: int = state_variable(unit='adim', description="The blade's shape in its class (Lseed)",
                                     value_comment="-1: none", kind=DESCRIPTOR, default=-1,
                                     scale=SCALES.Organ, dtype='int')
    shape_mature_length: float = state_variable(unit='m', description="The blade's final length, its shape's length",
                                                value_comment="the phytomer table's Ll at its creation, then "
                                                              "final_blade_length once given",
                                                kind=DESCRIPTOR, default=0., scale=SCALES.Organ)
    shape_max_width: float = state_variable(unit='m', description="The blade's maximal width, its shape's width",
                                            value_comment="the phytomer table's Lw at its creation, then "
                                                          "final_blade_width once given",
                                            kind=DESCRIPTOR, default=0., scale=SCALES.Organ)

    organ_area: float = state_variable(unit='m2', description="The organ's visible area (its visible elements')",
                                       kind=DESCRIPTOR, default=0., scale=SCALES.Organ)

    # --- The elements ---
    element_index: int = state_variable(unit='adim', description="The element's position in its organ, from 0",
                                        kind=DESCRIPTOR, default=0, scale=SCALES.SubOrgan,
                                        dtype='int')
    length: float = state_variable(unit='m', description="The element's length", kind=DESCRIPTOR,
                                   default=0., scale=SCALES.SubOrgan)
    area: float = state_variable(unit='m2', description="The element's area", kind=DESCRIPTOR,
                                 default=0., scale=SCALES.SubOrgan)
    green_length: float = state_variable(unit='m', description="The element's green length",
                                         kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    green_area: float = state_variable(unit='m2', description="The element's green area",
                                       kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    senesced_length: float = state_variable(unit='m', description="The element's senesced length",
                                            kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    senesced_area: float = state_variable(unit='m2', description="The element's senesced area",
                                          kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    is_green: int = state_variable(unit='adim', description="1 when the element is more green than senesced",
                                   kind=DESCRIPTOR, default=1, scale=SCALES.SubOrgan, dtype='int')
    rolled_part: float = state_variable(unit='m', description="The blade element's rolled length",
                                        kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    rolled_diameter: float = state_variable(unit='m', description="The blade element's rolled part's diameter",
                                            kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    width: float = state_variable(unit='m', description="The blade element's mean width",
                                  kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    flat_base: float = state_variable(unit='adim',
                                      description="The blade element's flat part's base, relative to the blade's "
                                                  "visible length (Adel's srb)", value_comment="NaN: no flat part",
                                      kind=DESCRIPTOR, default=np.nan, scale=SCALES.SubOrgan)
    flat_top: float = state_variable(unit='adim',
                                     description="The blade element's flat part's top, relative to the blade's visible "
                                                 "length (Adel's srt)", value_comment="NaN: no flat part",
                                     kind=DESCRIPTOR, default=np.nan, scale=SCALES.SubOrgan)

    # --- The elements' placement (the turtle's frame where each is drawn, and its end) ---
    x: float = planting_variable("x")
    y: float = planting_variable("y")
    z: float = planting_variable("z")
    rotation: float = planting_variable("rotation")
    x1: float = state_variable(unit='m', description="The element's base, x (topological: its parent's end)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    y1: float = state_variable(unit='m', description="The element's base, y (topological: its parent's end)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    z1: float = state_variable(unit='m', description="The element's base, z (topological: its parent's end)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    x2: float = state_variable(unit='m', description="The element's end, x (topological: its base plus its chord)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    y2: float = state_variable(unit='m', description="The element's end, y (topological: its base plus its chord)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    z2: float = state_variable(unit='m', description="The element's end, z (topological: its base plus its chord)",
                               kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    # where the turtle draws them (a blade element's along its blade; a hidden one on the visible parts below it)
    display_x1: float = display_end("x1")
    display_y1: float = display_end("y1")
    display_z1: float = display_end("z1")
    display_x2: float = display_end("x2")
    display_y2: float = display_end("y2")
    display_z2: float = display_end("z2")
    frame_x: float = state_variable(unit='m', description="The origin of the element's frame, x: where it is drawn "
                                    "from (a stem element's base; a blade element's, its rolled part's base)",
                                    kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    frame_y: float = state_variable(unit='m', description="The origin of the element's frame, y",
                                    kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    frame_z: float = state_variable(unit='m', description="The origin of the element's frame, z",
                                    kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    heading_x: float = state_variable(unit='adim', description="The element's heading, x",
                                      kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    heading_y: float = state_variable(unit='adim', description="The element's heading, y",
                                      kind=DESCRIPTOR, default=0., scale=SCALES.SubOrgan)
    heading_z: float = state_variable(unit='adim', description="The element's heading, z",
                                      kind=DESCRIPTOR, default=1., scale=SCALES.SubOrgan)
    up_x: float = state_variable(unit='adim',
                                 description="The element's up, x (its blade's side: the shape's x is along it)",
                                 kind=DESCRIPTOR, default=-1., scale=SCALES.SubOrgan)
    up_y: float = state_variable(unit='adim', description="The element's up, y", kind=DESCRIPTOR,
                                 default=0., scale=SCALES.SubOrgan)
    up_z: float = state_variable(unit='adim', description="The element's up, z", kind=DESCRIPTOR,
                                 default=0., scale=SCALES.SubOrgan)

    # --- Parameters ---
    axes_table: object = parameter(unit='',
                                   description="Adel's axeT: one row per axis of its plants (plant, axe, nf, HS_final, "
                                               "hasEar, ...)", value_comment="each plant's (its scenario's)", by=BY,
                                   default=None, scale=SCALES.Plant, dtype='object')
    phytomer_table: object = parameter(unit='',
                                       description="Adel's phytoT: one row per plant, axis and rank (plant, axe, n, "
                                                   "Ll, Lw, Gl, Gd, El, Ed, Azim, Lindex, Lseed)",
                                       value_comment="each plant's (its scenario's)", by=BY, default=None,
                                       scale=SCALES.Plant, dtype='object')
    table_plant: int = parameter(unit='adim', description="The plant's number in the tables (their column 'plant')",
                                 by=BY, default=1)
    n_sect: int = setting(unit='adim', description="The number of elements of a blade", default=1)
    length_unit: float = parameter(unit='m', description="The tables' length unit in m (Adel's cm: 0.01)", by=BY,
                                   default=0.01)
    leaves: object = setting(description="The leaf shape database (LeafShapes; None: Adel's samples)", default=None)
    min_length: float = setting(unit='m',
                                description="The shortest element drawn (Adel's min_length, 0.01 cm); a blade's "
                                            "triangles smaller than its square / 100 are left out", default=1e-4)
    tiller_inclination: int = parameter(unit='adim',
                                        description="1: the axes inclined by static Adel's rules (the axes table's "
                                                    "incT, dredT, azTb, azT; Adel.R's getIncT, axe_inclination); 0: "
                                                    "vertical, as AdelDyn", by=BY, default=0)
    own_axis_top: int = parameter(unit='adim',
                                  description="1: at a change of axis, the turtle starts on the top the element's "
                                              "topological parent left: a tiller of a tiller on its bearer's phytomer, "
                                              "an axis continuing after its own tiller on its own top (correction "
                                              "AD-37); 0: Adel's turtle, which continues it on the main stem's top",
                                  by=BY, default=0)
    inclination_start_haun_stage: float = parameter(unit='adim',
                                                    description="The Haun stage from which a tiller's inclination "
                                                                "ramps (Adel's HSstart_inclination_tiller)", by=BY,
                                                    default=1.)
    inclination_rate: float = parameter(unit='°',
                                        description="A tiller's inclination gained per Haun stage (Adel's "
                                                    "rate_inclination_tiller)", by=BY, default=30.)
    phyllochron: float = parameter(unit='°Cd',
                                   description="The main stem's phyllochron, for its blades' age classes in a dynamic "
                                               "leaf database", by=BY, default=88.92)
    tiller_phyllochron: float = parameter(unit='°Cd',
                                          description="The tillers' phyllochron, for their blades' age classes in a "
                                                      "dynamic leaf database (Adel's T1's for every tiller)", by=BY,
                                          default=91.6)
    leaf_number: float = parameter(unit='adim',
                                   description="Without tables: the main stem's number of leaves of a generated plant "
                                               "(tables.wheat_tables)",
                                   value_comment="NaN: 12 when another scalar is given; none given: Adel's sample "
                                                 "plant", by=BY, default=np.nan)
    tiller_number: float = parameter(unit='adim', description="Without tables: a generated plant's number of tillers",
                                     value_comment="NaN: 4", by=BY, default=np.nan)
    longest_blade_length: float = parameter(unit='m',
                                            description="Without tables: a generated plant's longest blade's length",
                                            value_comment="NaN: 0.19", by=BY, default=np.nan)
    longest_blade_width: float = parameter(unit='m',
                                           description="Without tables: a generated plant's longest blade's width",
                                           value_comment="NaN: 0.016", by=BY, default=np.nan)

    # --- The plant's structure ---

    @classmethod
    def initiate_plant(cls, g, plant, parameters):
        """The plant's tables (tables.plant_tables: its shoot_tables), its main stem (its Axis, GrowthUnit) and its
        first phytomer; its axes' rows are the population's tables' (its tables at their offset among the plants built
        before it).

        With a scenario's `shoot_anchor` (a selector of the vertex another initiator built, given to metafspm's
        MPG.find_anchor: {"label": ..., "scale": ..., "where": {...}, "first" or "last": ...}, e.g. a crown's top
        StemElement, {"label": "StemElement", "last": True}), the main stem's first element follows that vertex ('<')
        instead of starting the plant's topology; it is recorded as the plant's shoot_base (the elements it bears are
        the architecture's: SHOOT_ELEMENTS), and the anchor's z2, when its model gave it, as the plant's
        shoot_base_height (the main stem's base: the anchor's top)."""
        axes, phytomers, number = plant_tables(parameters, dict(table_plant=cls.table_plant,
                                                                phyllochron=cls.phyllochron,
                                                                length_unit=cls.length_unit))
        selector = parameters.get("shoot_anchor")
        anchored_on = None if selector is None else g.find_anchor(plant, **selector)
        record_plant_tables(g, plant, (axes, phytomers, number))
        built = g.properties()["shoot_tables"]
        offset = int(population_tables([built[p] for p in sorted(built) if p <= plant])[3][-1])
        n_sect = int(parameters.get("n_sect", cls.n_sect))
        unit = float(parameters.get("length_unit", cls.length_unit))          # the plant's
        rows = axes.reset_index(drop=True)
        mine = plant_rows(rows, number)
        main_stem = mine[0]
        scales, labels = g.scales, g.labels
        axis = g.add_component(plant, **PropsConfig(scale=scales.Axis, edge_type='/', label=labels.Axis.Shoot))
        unit_vid = g.add_component(axis, **PropsConfig(scale=scales.GrowthUnit, edge_type='/',
                                                       label=labels.GrowthUnit.Shoot))
        values = cls._axis_values(rows, main_stem, unit)
        values["axis_row"] += offset                              # the population's rows
        values["bearer_row"] = values["bearer_row"] + offset if values["bearer_row"] >= 0 else -1
        for name, value in dict(values, axis_phytomers=1).items():
            g.property(name)[axis] = value
        dimensions = phytomer_dimensions(phytomers, rows, 1.)                  # in the table's unit
        complexes, elements = phytomer_layout(labels, n_sect)
        created = []
        for anchor, scale, label in complexes:
            parent = unit_vid if anchor == "unit" else created[anchor]
            created.append(g.add_component(parent, **PropsConfig(scale=scale, edge_type='/', label=label)))
        cls._set_phytomer(g, created, dimensions, main_stem, 1, int(rows.loc[main_stem, "nf"]), unit)
        vids = []
        for complex_index, parent, edge, label, _, index in elements:
            if parent == "previous" and anchored_on is not None:          # on the anchor's top
                vid = g.add_component_with_topo(created[complex_index], anchored_on, **PropsConfig(
                    scale=scales.SubOrgan, edge_type='<', label=label))
                g.property("shoot_base")[plant] = vid
                top = g.properties().get("z2", {}).get(anchored_on)
                if top is not None:
                    g.property("shoot_base_height")[plant] = float(top)
            elif parent == "previous":
                vid = g.add_component(created[complex_index], **PropsConfig(scale=scales.SubOrgan, edge_type='/',
                                                                            label=label))
            else:
                vid = g.add_component_with_topo(created[complex_index], vids[parent], **PropsConfig(
                    scale=scales.SubOrgan, edge_type=edge, label=label))
            g.property("element_index")[vid] = index
            vids.append(vid)

    @staticmethod
    def _axis_values(rows, row, unit):
        """An axis' values from its row of the axes table (its bearer's row among the same plant's); lengths in m
        (*unit*: the table's unit)."""
        name = rows.loc[row, "axe"]
        bearer, rank = bearer_of(name)
        same_plant = rows.index[(rows["plant"] == rows.loc[row, "plant"]) & (rows["axe"] == bearer)].tolist()
        haun_stage = rows.loc[row, "HS_final"]

        def angle(column):                             # optional columns: 0 when missing
            value = rows.loc[row, column] if column in rows.columns else np.nan
            return float(value) if pd.notna(value) else 0.

        return dict(axis_row=row, bearer_row=same_plant[0] if same_plant else -1, insertion_rank=rank,
                    final_leaf_number=int(rows.loc[row, "nf"]), has_ear=int(rows.loc[row, "hasEar"]),
                    final_haun_stage=float(haun_stage) if pd.notna(haun_stage) else float(rows.loc[row, "nf"]),
                    table_inclination=angle("incT"), straightening_distance=angle("dredT") * unit,
                    insertion_azimuth=angle("azTb"), first_azimuth=angle("azT"))

    @classmethod
    def _set_phytomer(cls, g, created, dimensions, row, rank, nff, unit):
        """A phytomer's and its organs' values on the MPG (initiate_plant's)."""
        for name, value in cls._phytomer_values(dimensions, np.array([row]), np.array([rank]), np.array([nff]),
                                                np.array([unit])).items():
            level, values = value
            for vid, v in zip([created[i] for i in level], values):
                g.property(name)[vid] = np.asarray(v).reshape(-1)[0].item()     # this one phytomer's

    @staticmethod
    def _phytomer_values(dimensions, rows, ranks, nff, units):
        """{variable: (complex positions in a phytomer's layout, values per phytomer per position)}: for phytomers
        of axis rows *rows* at ranks *ranks*, from *dimensions* in the tables' units (*units*: each phytomer's plant's
        length_unit, m); layout positions 0 phytomer, 1 internode, 2 sheath, 3 blade."""
        def at(name):
            return dimensions[name][rows, ranks]

        def length(name):
            return at(name) * units
        shape_class, shape_seed = at("Lindex"), at("Lseed")
        valid = (shape_class != -999) & (shape_seed != -999) & np.isfinite(shape_class) & np.isfinite(shape_seed)
        ntop = nff - ranks + 1
        return {
            "phytomer_rank": ([0], [ranks]),
            "organ_rank": ([1, 2, 3], [ranks, ranks, ranks]),
            "ntop": ([1, 2, 3], [ntop, ntop, ntop]),
            "diameter": ([1, 2], [length("Ed"), length("Gd")]),
            "organ_azimuth": ([1], [at("Azim")]),
            "shape_mature_length": ([3], [length("Ll")]),
            "shape_max_width": ([3], [length("Lw")]),
            "shape_class": ([3], [np.where(valid, shape_class, -1).astype(int)]),
            "shape_seed": ([3], [np.where(valid, shape_seed, -1).astype(int)]),
        }

    def tables(self):
        """(axes table, phytomer table, per plant its number in them, per plant its rows' offset): the population's
        tables (tables.population_tables), the ones axis_row indexes."""
        return population_tables(self.data_structure.get("shoot_tables"))

    def axis_names(self):
        """Each Axis' name (MS, T1, ...), from its row in the population's axes table ("" for another model's axis,
        without a row: a crown's)."""
        rows = np.asarray(self.data_structure.get("axis_row"))
        names = self.tables()[0]["axe"].to_numpy()[np.maximum(rows, 0)]
        return np.where(rows >= 0, names, "")

    def organ_names(self):
        """Each Organ's kind: "internode", "sheath" or "blade"."""
        _, organ_kind = self._kinds()
        ds = self.data_structure
        kind = np.full(len(ds.entity_ids("Organ")), -1)
        kind[np.asarray(ds.owner("Organ"))] = organ_kind
        return np.array(ORGANS + ("",), dtype=object)[kind]

    def _dimension_arrays(self):
        """The population's phytomer table as arrays (tables.phytomer_dimensions), in the tables' units, read once."""
        def build():
            axes, phytomers = self.tables()[:2]
            return phytomer_dimensions(phytomers, axes, 1.)
        return self.cached("dimensions", build, key=id(self.tables()[0]))

    def _axes_lookup(self):
        """(rows of each table plant, the values of a row in a plant's length unit) of the axes table
        (tables.plant_rows, _axis_values), each read once."""
        def build():
            rows = self.tables()[0]
            return rows, {int(number): plant_rows(rows, int(number)) for number in set(rows["plant"])}, {}
        rows, rows_of, values = self.cached("axes", build, key=id(self.tables()[0]))

        def values_of(row, unit):
            if (row, unit) not in values:
                values[(row, unit)] = self._axis_values(rows, row, unit)
            return values[(row, unit)]
        return rows_of, values_of

    @segmentation
    def _axis_phytomers(self, axis_phytomers, phytomer_number, axis_row) -> dict:
        """The axes the plants' axis_number asks for (in their axes table's order, each with its first phytomer, once
        its bearer has the phytomer it is inserted on), and the phytomers the axes' phytomer_number asks for: created
        with their organs and elements (NewNodes), their dimensions and shapes from the phytomer table. A call creates
        all of them, in order: the new axes, the phytomers of the axes there before (as many per axis as asked, each
        after the previous one), then the axes waiting for a bearer's phytomer this call created. A new axis has its
        first phytomer: the development gives its phytomer_number after the call's segmentation
        (ScenarioShoot.grow then segments again). Its axes are those with a row in the axes table (another
        model's, a crown's, has none: left as they are)."""
        ds = self.data_structure
        g = ds.structure
        plant_rows, axis_values = self._axes_lookup()
        plant_unit = np.asarray(ds.parameter_view("length_unit", "Plant"), dtype=float)        # each plant's
        count = np.array(axis_phytomers, dtype=int)
        target = np.floor(np.asarray(phytomer_number, dtype=float)).astype(int)
        axes = np.asarray(ds.entity_ids("Axis"))
        axis_of = {(g.complex(int(axis)), int(row)): int(axis) for axis, row in zip(axes, axis_row) if row >= 0}
        batch = _Batch(g.labels, int(self.n_sect))
        plants = np.asarray(ds.entity_ids("Plant"))
        wanted = np.floor(np.asarray(ds.get("axis_number"), dtype=float)).astype(int)
        table_plant = self.tables()[2]

        def new_axes():
            """The axes asked for, in order, with their first phytomer, up to the first whose bearer's phytomer is
            not there yet (in the MPG or the batch)."""
            for plant, number, table, unit in zip(plants.tolist(), wanted.tolist(), table_plant.tolist(),
                                                  plant_unit.tolist()):
                mine = plant_rows[table]
                existing = [row for row in mine if (plant, row) in axis_of or (plant, row) in batch.internodes]
                for row in mine[len(existing):number]:
                    values = axis_values(row, unit)
                    predecessor = self._predecessor(g, None, 1, axis_of, plant, values["bearer_row"],
                                                    values["insertion_rank"], batch)
                    if predecessor is None:
                        break                               # its bearer's phytomer is not there yet
                    batch.axis(plant, row, dict(values, axis_phytomers=1))
                    batch.phytomer(None, predecessor, plant, row, 1, values["final_leaf_number"], unit,
                                   branched=True)

        new_axes()

        # New phytomers on the axes
        axis_plant = ds.index_of(np.array([g.complex(axis) for axis in axes.tolist()], dtype=np.int64), "Plant")
        for a, axis in enumerate(axes.tolist()):
            if axis_row[a] < 0:                             # another model's axis
                continue
            unit = float(plant_unit[axis_plant[a]])
            values = axis_values(int(axis_row[a]), unit)
            for rank in range(count[a] + 1, target[a] + 1):
                predecessor = self._predecessor(g, axis, rank, axis_of, g.complex(axis), values["bearer_row"],
                                                values["insertion_rank"], batch, row=int(axis_row[a]))
                batch.phytomer(g.components(axis)[0], predecessor, g.complex(axis), int(axis_row[a]), rank,
                               values["final_leaf_number"], unit, branched=False)
                count[a] = rank

        new_axes()                                          # those on a bearer's phytomer just created
        return {"axis_phytomers": count,
                "new_phytomers": batch.new_nodes(self._dimension_arrays(), self.variable_specs)}

    # --- The elements ---

    def _kinds(self):
        """Per node: its element kind (0 hidden, 1 stem, 2 leaf) and its organ's (0 internode, 1 sheath, 2 blade),
        from the MPG's labels (ds.labels), once per topology."""
        def build():
            ds = self.data_structure
            element = {"HiddenElement": 0, "StemElement": 1, "LeafElement": 2}
            organ = {"StemInternode": 0, "Sheath": 1, "Blade": 2}
            return (np.array([element.get(label, -1) for label in ds.labels().tolist()]),
                    np.array([organ.get(label, -1) for label in ds.labels(location="Organ").tolist()]))
        return self.cached("kinds", build)

    def _leaf_shapes(self):
        return default_leaves() if self.leaves is None else self.leaves

    @postsegmentation
    def _elements(self, element_index, organ_length, visible_length, organ_senesced_length, rolled_length, diameter,
                  shape_mature_length, shape_max_width, final_blade_length, final_blade_width, shape_class,
                  organ_area) -> dict:
        """The elements from their organs' lengths (Adel's internode_elements, sheath_elements, blade_elements):
        a stem organ's hidden element (its length out of sight) and visible one, its green and senesced parts; a
        blade's hidden part, and its sectors along its mature shape (its visible part the shape's top; its final
        dimensions the development's once given), each one's green and senesced parts, rolled part, flat part and
        width; areas from the diameters and the shapes' width profiles; the organs' visible areas (another model's
        organs, without its elements, keep theirs)."""
        shape_mature_length = np.where(np.isfinite(final_blade_length), final_blade_length, shape_mature_length)
        shape_max_width = np.where(np.isfinite(final_blade_width), final_blade_width, shape_max_width)
        ds = self.data_structure
        element_kind, organ_kind = self._kinds()
        organ = np.asarray(ds.owner("Organ"))
        n = element_kind.size
        l, l_visible = organ_length[organ], visible_length[organ]
        l_senesced, l_rolled = organ_senesced_length[organ], rolled_length[organ]
        out = {name: np.zeros(n) for name in ("area", "green_length", "green_area", "senesced_length",
                                              "senesced_area", "rolled_part", "rolled_diameter", "width")}
        out["length"] = np.zeros(n)
        out["flat_base"], out["flat_top"] = np.full(n, np.nan), np.full(n, np.nan)
        is_green = np.ones(n, dtype=int)
        hidden = np.maximum(l - l_visible, 0.)

        # Stem organs (internodes, sheaths): a cylinder of the organ's diameter
        stem = organ_kind < 2
        d = diameter[organ]
        green = l_visible - np.minimum(l_senesced, l_visible)
        senesced = l_visible - green
        these = stem & (element_kind == 0)
        out["length"][these] = hidden[these]
        out["area"][these] = np.pi * hidden[these] * d[these]
        these = stem & (element_kind == 1)
        out["length"][these] = l_visible[these]
        out["area"][these] = np.pi * l_visible[these] * d[these]
        out["green_length"][these], out["senesced_length"][these] = green[these], senesced[these]
        out["green_area"][these] = np.pi * green[these] * d[these]
        out["senesced_area"][these] = np.pi * senesced[these] * d[these]
        is_green[these] = green[these] >= senesced[these]

        # Blades: the visible part is the mature shape's top
        shapes = self._leaf_shapes()
        blade = organ_kind == 2
        mature, max_width = shape_mature_length[organ], shape_max_width[organ]
        classes = shape_class[organ]
        has_shape = blade & (classes >= 1) & (mature > 0)
        these = blade & (element_kind == 0)
        out["length"][these] = hidden[these]
        shaped = these & has_shape & (hidden > 0)
        out["area"][shaped] = shapes.relative_area(classes[shaped], (mature[shaped] - l[shaped]) / mature[shaped],
                                                   (mature[shaped] - l_visible[shaped]) / mature[shaped]) * mature[
            shaped] * max_width[shaped]
        sector = blade & (element_kind == 2) & (mature > 0) & (l_visible > 1e-6)     # without a shape: no area
        if sector.any():
            self._blade_sectors(out, is_green, sector, element_index, l_visible, l_senesced, l_rolled, mature,
                                max_width, classes, organ, shapes)

        # The organs' visible areas
        visible = ((element_kind == 1) & stem) | ((element_kind == 2) & blade)
        n_organs = len(ds.entity_ids("Organ"))
        areas = np.bincount(organ[visible], weights=out["area"][visible], minlength=n_organs)
        areas = np.where(np.bincount(organ, minlength=n_organs) > 0, areas, organ_area)
        return dict(out, is_green=is_green, organ_area=areas, shape_mature_length=shape_mature_length,
                    shape_max_width=shape_max_width)

    @staticmethod
    def _blade_sectors(out, is_green, sector, element_index, l_visible, l_senesced, l_rolled, mature, max_width,
                       classes, organ, shapes):
        """Each visible blade's sectors (Adel's blade_elements): sector k covers [(k - 1), k] · L / n_sect of the
        mature shape (L its length), its visible part above L - visible length, its senesced part above L -
        senesced length, its flat part above L - (visible - rolled) length. A blade without a shape class has its
        sectors' lengths and flat parts, no area and no width (Adel's, its shape key None)."""
        k = element_index[sector]
        n_sect = np.bincount(organ[sector], minlength=organ.max() + 1)[organ[sector]]
        L, w = mature[sector], max_width[sector]
        visible = l_visible[sector]
        rolled = np.clip(l_rolled[sector], 0., visible)
        rolled = np.where(rolled < 1e-6, 0., rolled)
        flat = visible - rolled
        green = visible - np.minimum(l_senesced[sector], visible)
        senesced = visible - green
        limit_visible, limit_senesced, limit_rolled = L - visible, L - senesced, L - flat
        step = L / n_sect
        top = k * step
        visible_part = np.minimum(step, np.maximum(0., top - limit_visible))
        base = top - visible_part
        green_top = np.minimum(top, np.maximum(base, limit_senesced))
        green_part, senesced_part = green_top - base, top - green_top
        shaped = classes[sector] >= 1
        c = np.where(shaped, classes[sector], 1)
        green_area = np.where(shaped & (green_part > 0), shapes.relative_area(c, base / L, green_top / L) * L * w, 0.)
        senesced_area = np.where(shaped & (senesced_part > 0), shapes.relative_area(c, green_top / L, top / L) * L * w,
                                 0.)
        total = green_area + senesced_area
        flat_part = np.minimum(visible_part, np.maximum(0., top - limit_rolled))
        has_flat = flat_part > 0
        flat_bottom = top - flat_part
        flat_green_top = np.minimum(top, np.maximum(flat_bottom, limit_senesced))
        with np.errstate(divide="ignore", invalid="ignore"):
            rolled_part = np.where(flat_part < visible_part, visible_part - flat_part, 0.)
            rolled_diameter = np.where(rolled_part > 0, total * rolled_part / visible_part / np.pi / rolled_part, 0.)
            green_width = np.where(shaped & (green_part > 0), shapes.mean_width(c, base / L, green_top / L) * w, 0.)
            senesced_width = np.where(shaped & (senesced_part > 0), shapes.mean_width(c, green_top / L, top / L) * w,
                                      0.)
            mean_width = np.maximum(green_width, senesced_width)
        index = np.flatnonzero(sector)
        visible_any = visible_part > 0
        out["length"][index] = np.where(visible_any, senesced_part + green_part, 0.)
        out["area"][index] = np.where(visible_any, total, 0.)
        out["green_length"][index], out["green_area"][index] = green_part, green_area
        out["senesced_length"][index], out["senesced_area"][index] = senesced_part, senesced_area
        out["rolled_part"][index], out["rolled_diameter"][index] = rolled_part, rolled_diameter
        out["width"][index] = mean_width
        out["flat_base"][index] = np.where(has_flat, (flat_bottom - limit_visible) / visible, np.nan)
        out["flat_top"][index] = np.where(has_flat, (top - limit_visible) / visible, np.nan)
        is_green[index] = green_part > senesced_part

    # --- The placement and the triangles ---

    def _traversal(self):
        """Per node: its position in Adel's traversal and its frame's source (turtle.traversal_order), whether it is
        its plant's base (without a topological parent among the elements: none, or the vertex it is anchored on),
        and its topological parent's index (-1: none among them); once per topology."""
        def build():
            ds = self.data_structure
            g = ds.structure
            nodes = ds.entity_ids("node").tolist()
            index = {vid: i for i, vid in enumerate(nodes)}
            parents = [g.parent(v) for v in nodes]
            return (*traversal_order(g, nodes), np.array([parent not in index for parent in parents]),
                    np.array([index.get(parent, -1) for parent in parents], dtype=np.int64))
        return self.cached("traversal", build)

    @postsegmentation
    def _first_inclination(self, first_inclination, axis_azimuth, table_inclination, bearer_row, haun_stage,
                           insertion_azimuth, axis_row, tiller_inclination, inclination_start_haun_stage,
                           inclination_rate) -> dict:
        """With tiller_inclination, each axis' first internode's inclination by static Adel's rule
        (turtle.base_inclination: a tiller's ramps with its Haun stage) and its azimuth at its insertion, its
        insertion azimuth (azTb); an axis' step, its parameters its plant's. Without it, nothing is written; nor on
        another model's axis (without a row in the axes table)."""
        on = np.asarray(tiller_inclination) > 0
        if not on.any():
            return None
        on = on & (np.asarray(axis_row) >= 0)
        first = base_inclination(np.asarray(table_inclination, dtype=float), np.asarray(bearer_row) < 0,
                                 np.asarray(haun_stage, dtype=float),
                                 np.asarray(inclination_start_haun_stage, dtype=float),
                                 np.asarray(inclination_rate, dtype=float))
        return {"first_inclination": np.where(on, first, first_inclination),
                "axis_azimuth": np.where(on, np.asarray(insertion_azimuth, dtype=float), axis_azimuth)}

    @postsegmentation
    def _inclination(self, inclination, organ_azimuth, visible_length, organ_rank, first_inclination,
                     straightening_distance, first_azimuth, tiller_inclination) -> dict:
        """With tiller_inclination, the stem organs' inclinations by static Adel's rule (kernels.stem_inclinations):
        each axis' first internode inclined by its first inclination, a tiller straightened at its straightening
        distance; its first internode's azimuth (azT). An organ's step, its parameter its plant's. Without it,
        nothing is written (AdelDyn's: 0, the phytomer table's Azim)."""
        on = np.asarray(tiller_inclination) > 0
        if not on.any():
            return None
        ds = self.data_structure
        _, organ_kind = self._kinds()
        kind = np.full(len(inclination), -1)
        kind[np.asarray(ds.owner("Organ"))] = organ_kind
        axis = ds.owner("Axis", of="Organ")
        stem = np.flatnonzero((kind == 0) | (kind == 1))
        stem = stem[np.lexsort((kind[stem], organ_rank[stem], axis[stem]))]
        inclined = np.array(inclination, dtype=float)
        computed = stem_inclinations(axis[stem], kind[stem].astype(np.int64),
                                     np.asarray(visible_length, dtype=float)[stem],
                                     np.asarray(first_inclination, dtype=float),
                                     np.asarray(straightening_distance, dtype=float))
        chosen = on[stem]
        inclined[stem[chosen]] = computed[chosen]
        azimuth = np.array(organ_azimuth, dtype=float)
        first_internode = (kind == 0) & (np.asarray(organ_rank) == 1) & on
        azimuth[first_internode] = np.asarray(first_azimuth, dtype=float)[axis[first_internode]]
        return {"inclination": inclined, "organ_azimuth": azimuth}

    @postsegmentation
    def _placement(self, length, rolled_part, inclination, organ_azimuth, axis_azimuth, bearer_row,
                   phytomer_rank, x, y, z, rotation, shoot_base_height, own_axis_top) -> dict:
        """The elements' placement, their frames by Adel's turtle (kernels.interpret): each one's frame (its origin,
        heading and up) where it is drawn, its base and its end; from its plant's position and rotation (its main
        stem's base at its shoot_base_height when given: an anchored main stem's, on the anchor's top), on the
        elements' lengths. A stem element's base is its frame's origin, its end that plus its length along its
        heading (a hidden element's: its base). A blade element is drawn in its blade's frame (Adel's: a blade's
        elements from its collar, each moved up by the rolled parts below it): its base is its rolled part's base
        when it has one, else its flat part's midrib's first point; its end its flat part's midrib's last point, its
        rolled part's end without a flat part drawn (kernels.blade_element_midrib_ends): its display ends. Its
        topological ends (x1 ... z2) chain the elements' chords from the plant's base, each element's base at its
        parent's end: a stem or hidden element's chord its heading times its length, a blade element's its drawn
        chord (its drawn end minus its drawn base)."""
        ds = self.data_structure
        element_kind, organ_kind = self._kinds()
        order, frame_source, plant_base, parent = self._traversal()
        organ, axis = np.asarray(ds.owner("Organ")), np.asarray(ds.owner("Axis"))
        plant, phytomer = np.asarray(ds.owner("Plant")), np.asarray(ds.owner("Phytomer"))
        plant_position = np.column_stack([np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                          np.where(np.isnan(shoot_base_height), np.asarray(z, dtype=float),
                                                   shoot_base_height)]).reshape(-1, 3)    # the planting's
        plant_rotation = np.asarray(rotation, dtype=float)
        axis_is_main = np.asarray(bearer_row) < 0
        stem = element_kind == 1
        position, heading, up, end = interpret(
            order, element_kind.astype(np.int64), plant_base, axis.astype(np.int64), axis_is_main,
            (np.asarray(phytomer_rank)[phytomer] == 1) & ~axis_is_main[axis], plant_position[plant],
            plant_rotation[plant], np.asarray(axis_azimuth, dtype=float)[axis],
            np.where(stem, np.asarray(inclination, dtype=float)[organ], 0.),
            np.where(stem & (organ_kind == 0), np.asarray(organ_azimuth, dtype=float)[organ], 0.),
            np.asarray(length, dtype=float), np.asarray(rolled_part, dtype=float), frame_source, parent,
            np.asarray(own_axis_top) > 0)
        base = position.copy()
        nodes, ends = self._blade_midrib_ends(element_kind, organ, np.asarray(length), np.asarray(rolled_part))
        if nodes.size:
            placed = place(ends[:, None], position[nodes], heading[nodes], up[nodes])[:, 0]
            rolled = np.asarray(rolled_part, dtype=float)[nodes] > 0
            base[nodes] = np.where(rolled[:, None], position[nodes], placed[:, 0])
            end[nodes] = placed[:, 1]
        chord = np.where((element_kind == 2)[:, None], end - base, heading * np.asarray(length, dtype=float)[:, None])
        root = plant_position[plant]                          # each element's plant's base, for its first element
        ds = self.data_structure
        top = np.column_stack([ds.accumulate(chord[:, axis], direction="down") for axis in range(3)]) + root
        bottom = top - chord
        return {"x1": bottom[:, 0], "y1": bottom[:, 1], "z1": bottom[:, 2], "x2": top[:, 0], "y2": top[:, 1],
                "z2": top[:, 2], "display_x1": base[:, 0], "display_y1": base[:, 1], "display_z1": base[:, 2],
                "display_x2": end[:, 0], "display_y2": end[:, 1], "display_z2": end[:, 2],
                "frame_x": position[:, 0], "frame_y": position[:, 1], "frame_z": position[:, 2],
                "heading_x": heading[:, 0], "heading_y": heading[:, 1], "heading_z": heading[:, 2],
                "up_x": up[:, 0], "up_y": up[:, 1], "up_z": up[:, 2]}

    def frames(self):
        """(origin, heading, up), each (n, 3): the elements' frames where they are drawn."""
        ds = self.data_structure
        names = (("frame_x", "frame_y", "frame_z"), ("heading_x", "heading_y", "heading_z"), ("up_x", "up_y", "up_z"))
        return tuple(np.stack([np.asarray(ds.get(name), dtype=float) for name in frame], axis=1) for frame in names)

    def triangles(self, where=None) -> Triangles:
        """The elements' triangles (metafspm's coupling.geometry.Triangles), as Adel draws them: a stem element's
        slim cylinder of its organ's diameter; a blade element's rolled part's, of its rolled diameter, then its
        flat part's mesh above it (kernels' blade meshes). *where*: a node mask. Each triangle carries its "green"
        attribute: the share of its area in its element's green part, the element's first green_length from its base
        (the senescence going down from the top: a blade's sectors, a sheath's or an internode's visible part), its
        area beyond it the senesced part ("senesced": 1 - green)."""
        ds = self.data_structure
        element_kind, _ = self._kinds()
        organ = np.asarray(ds.owner("Organ"))
        length, rolled = np.asarray(ds.get("length")), np.asarray(ds.get("rolled_part"))
        minimum = self.min_length
        selected = np.ones(length.size, dtype=bool) if where is None else np.asarray(where, dtype=bool)
        leaf = selected & (element_kind == 2) & (length > 0)
        stem = selected & (element_kind == 1) & (length > 0) & (length >= minimum)
        cylinder = np.flatnonzero(stem | (leaf & (rolled > 0)))
        cylinder_length = np.where(stem, length, rolled)[cylinder]
        radius = np.where(stem, 0.5 * np.asarray(ds.get("diameter"))[organ],
                          0.5 * np.asarray(ds.get("rolled_diameter")))[cylinder]
        position, heading, up = self.frames()
        local = slim_cylinders(cylinder_length, radius, radius)
        vertices = place(local, position[cylinder], heading[cylinder], up[cylinder]).reshape(-1, 3, 3)
        owner = np.repeat(cylinder, 8)
        green_length = np.asarray(ds.get("green_length"), dtype=float)
        green = self._green_shares(local.reshape(-1, 3, 3)[..., 2], green_length[owner], length[owner])
        blade_vertices, blade_owner, blade_abscissae = self._blade_triangles(leaf, rolled, minimum, organ)
        if blade_owner.size:
            placed = place(blade_vertices[:, None], position[blade_owner], heading[blade_owner],
                           up[blade_owner])[:, 0]
            vertices, owner = np.concatenate([vertices, placed]), np.concatenate([owner, blade_owner])
            green = np.concatenate([green, self._green_shares(blade_abscissae, green_length[blade_owner],
                                                              length[blade_owner])])
        return Triangles(vertices, owner, attributes={"green": green})

    @staticmethod
    def _green_shares(along, green_length, length) -> np.ndarray:
        """The share of each triangle's area within its element's first *green_length* of its *length*: its green
        part, from its vertices' distances *along* its element from its base ((T, 3): a cylinder's local z, a blade
        mesh's abscissa along its midrib above its rolled part), linear over the triangle (metafspm's
        Triangles.fraction_in, exact); a wholly green element's triangles 1 (its top cap included), a wholly senesced
        one's 0 (its base cap included)."""
        green_length = np.asarray(green_length, dtype=float)
        along = np.asarray(along, dtype=float)
        shifted = np.zeros(along.shape + (3,))
        shifted[..., 2] = along - green_length[:, None]
        shares = Triangles(shifted, np.zeros(len(shifted), dtype=np.int64)).fraction_in(HalfSpace(0., above=False))
        return np.where(green_length >= np.asarray(length, dtype=float), 1., np.where(green_length > 0., shares, 0.))

    def _blade_triangles(self, leaf, rolled, minimum, organ):
        """The flat parts' meshes of the blade elements *leaf* (a node mask) in their local frames, above their rolled
        parts: (T, 3, 3) vertices, (T,) node indices, (T, 3) their vertices' distances along their element from its
        base (the rolled part, then the midrib's abscissa). Each blade's shape is its class', seed's and, in a dynamic
        database, its age class' (its age in its axis' phyllochrons, minus 0.3, as Adel's). An element's mesh is
        kept from one call to the next while what it is made from is the same (its shape, visible length, flat
        part, ...): only the grown and new ones are meshed again, every one placed again."""
        ds = self.data_structure
        shapes = self._leaf_shapes()
        min_area = minimum ** 2 / 100.
        vertex_ids = ds.entity_ids("node")
        previous = self.kept("blade_meshes", {})
        cache = {}
        parts, owners, abscissae = [], [], []
        for node, made_from in self._meshed_blades(leaf, rolled, minimum, organ):
            made_from = made_from + (min_area,)
            vertex = int(vertex_ids[node])
            kept = previous.get(vertex)
            if kept is not None and kept[0] == made_from:
                triangles, along = kept[1], kept[2]
            else:
                triangles, along = blade_element_triangles(*shapes.fitted[made_from[1]], *made_from[2:8], min_area)
                triangles[..., 2] += rolled[node]
                along += rolled[node]
            cache[vertex] = (made_from, triangles, along)
            parts.append(triangles)
            owners.append(np.full(len(triangles), node, dtype=np.int64))
            abscissae.append(along)
        self.keep("blade_meshes", cache, saved=False)       # rebuilt when missing: not checkpointed
        if not parts:
            return np.zeros((0, 3, 3)), np.zeros(0, dtype=np.int64), np.zeros((0, 3))
        return np.concatenate(parts), np.concatenate(owners), np.concatenate(abscissae)

    def _meshed_blades(self, leaf, rolled, minimum, organ):
        """[(node, what its mesh is made from)]: the blade elements of *leaf* (a node mask) drawn with a mesh (a
        visible length of at least *minimum*, a flat part, a shape); made from (the shapes' id, its shape index, its
        inclination, its blade's mature length and maximal width, visible length, its flat part's base and top, its
        rolled part)."""
        ds = self.data_structure
        shapes = self._leaf_shapes()
        visible = np.asarray(ds.get("visible_length"))[organ]
        flat_base, flat_top = np.asarray(ds.get("flat_base")), np.asarray(ds.get("flat_top"))
        shaped = np.asarray(ds.get("shape_class"))[organ] >= 1
        meshed = np.flatnonzero(leaf & (visible >= minimum) & np.isfinite(flat_base) & shaped)
        if meshed.size == 0 or not shapes.fitted:
            return []
        blade = organ[meshed]
        main = np.asarray(ds.get("bearer_row"))[ds.owner("Axis", of="Organ")[blade]] < 0
        plant = np.asarray(ds.owner("Plant"))[meshed]
        phyllochron = np.where(main, np.asarray(ds.parameter_view("phyllochron", "Plant"))[plant],     # its plant's
                               np.asarray(ds.parameter_view("tiller_phyllochron", "Plant"))[plant])
        age_class = shapes.age_class(np.asarray(ds.get("age"))[blade] / phyllochron - 0.3)
        shape = shapes.shape_index(np.asarray(ds.get("shape_class"))[blade], np.asarray(ds.get("shape_seed"))[blade],
                                   age_class)
        inclination = 1. if shapes.dynamic else None          # a dynamic database's shapes carry their inclination
        mature, max_width = np.asarray(ds.get("shape_mature_length")), np.asarray(ds.get("shape_max_width"))
        blade_inclination = np.asarray(ds.get("inclination"))
        return [(node, (id(shapes), shape_index,
                        blade_inclination[organ_index] if inclination is None else inclination, mature[organ_index],
                        max_width[organ_index], visible[node], flat_base[node], flat_top[node], rolled[node]))
                for node, organ_index, shape_index in zip(meshed.tolist(), blade.tolist(), shape.tolist())
                if shape_index >= 0]

    def _blade_midrib_ends(self, element_kind, organ, length, rolled):
        """(nodes, ends (n, 2, 3)): the blade elements drawn with a mesh (triangles' leaf elements), their flat
        parts' midribs' first and last points in their frames, above their rolled parts; each kept from one call to
        the next while what its mesh is made from is the same."""
        shapes = self._leaf_shapes()
        vertex_ids = self.data_structure.entity_ids("node")
        previous = self.kept("blade_midrib_ends", {})
        cache = {}
        nodes, ends = [], []
        for node, made_from in self._meshed_blades((element_kind == 2) & (length > 0), rolled, self.min_length,
                                                   organ):
            vertex = int(vertex_ids[node])
            kept = previous.get(vertex)
            if kept is not None and kept[0] == made_from:
                points = kept[1]
            else:
                points = blade_element_midrib_ends(*shapes.fitted[made_from[1]], *made_from[2:8])
                points[:, 2] += rolled[node]
            cache[vertex] = (made_from, points)
            if len(points):
                nodes.append(node)
                ends.append(points)
        self.keep("blade_midrib_ends", cache, saved=False)  # rebuilt when missing: not checkpointed
        if not nodes:
            return np.zeros(0, dtype=np.int64), np.zeros((0, 2, 3))
        return np.asarray(nodes, dtype=np.int64), np.stack(ends)

    @staticmethod
    def _predecessor(g, axis, rank, axis_of, plant, bearer_row, insertion_rank, batch, row=None):
        """The element a new phytomer's internode follows: the previous internode's last element; for an axis' first
        phytomer, its bearer's phytomer's internode's first element (NO_PARENT for a main stem); None while that
        phytomer does not exist. An element of the MPG (its id), or of *batch* (a phytomer it holds: -1 - its index);
        *row*: the axis' row (for a phytomer after one in the batch)."""
        def internode_element(on_axis, on_row, on_rank, last):
            if on_axis is not None:
                phytomers = g.components(g.components(on_axis)[0])
                if len(phytomers) >= on_rank:
                    return g.components(g.components(phytomers[on_rank - 1])[0])[-1 if last else 0]
            return batch.internode_element(plant, on_row, on_rank, last)

        if rank > 1:
            return internode_element(axis, row, rank - 1, last=True)
        if bearer_row < 0:
            return NewNodes.NO_PARENT
        return internode_element(axis_of.get((plant, bearer_row)), bearer_row, max(insertion_rank, 1), last=False)


class _Batch:
    """The complexes and elements of a step's new axes and phytomers, gathered for one NewNodes."""

    def __init__(self, labels, n_sect):
        self.labels = labels
        self.complexes_layout, self.elements_layout = phytomer_layout(labels, n_sect)
        self.anchor, self.scale, self.label, self.complex_values = [], [], [], []
        self.parent, self.edge, self.element_label, self.complex_of, self.element_index = [], [], [], [], []
        self.phytomers = []                                 # (layout's first complex, axis row, rank, nff, unit)
        self.internodes = {}                    # {(plant, axis row): {rank: its internode's (first, last) elements}}
        self.last_unit = None

    def _complex(self, anchor, scale, label, values=None):
        self.anchor.append(anchor)
        self.scale.append(scale)
        self.label.append(label)
        self.complex_values.append(values or {})
        return len(self.anchor) - 1

    def axis(self, plant, row, values):
        """A new Axis of *plant* (its row *row*) with its GrowthUnit; the next phytomer is its first."""
        axis = self._complex(plant, SCALES.Axis, self.labels.Axis.Shoot, values)
        self.last_unit = self._complex(-1 - axis, SCALES.GrowthUnit, self.labels.GrowthUnit.Shoot)

    def internode_element(self, plant, row, rank, last):
        """The first (or *last*) element of the internode of the phytomer *rank* of *plant*'s axis *row* in the
        batch (-1 - its index), None when the batch does not hold it."""
        elements = self.internodes.get((plant, row), {}).get(rank)
        return None if elements is None else -1 - elements[1 if last else 0]

    def phytomer(self, unit, predecessor, plant, row, rank, nff, length_unit, branched):
        """A phytomer of *plant*'s axis *row* under the GrowthUnit *unit* (None: the new axis' just added), its
        internode after *predecessor* (an element's id, a batch element's -1 - its index, or NO_PARENT);
        *length_unit*: its plant's (m)."""
        unit = -1 - self.last_unit if unit is None else unit
        base = len(self.anchor)
        for where, scale, label in self.complexes_layout:
            self._complex(unit if where == "unit" else -1 - (base + where), scale, label)
        first = len(self.parent)
        for complex_index, where, link, label, _, index in self.elements_layout:
            if where == "previous":
                self.parent.append(predecessor)
                self.edge.append("+" if branched else "<")
            else:
                self.parent.append(-1 - (first + where))
                self.edge.append(link)
            self.element_label.append(label)
            self.complex_of.append(-1 - (base + complex_index))
            self.element_index.append(index)
        internode = [first + k for k, element in enumerate(self.elements_layout) if element[0] == 1]
        self.internodes.setdefault((plant, row), {})[rank] = (internode[0], internode[-1])
        self.phytomers.append((base, row, rank, nff, length_unit))

    def new_nodes(self, dimensions, specs) -> NewNodes:
        if not self.parent:
            return NewNodes(parent=np.zeros(0, dtype=np.int64))
        n = len(self.anchor)
        values = {}
        names = {name for complex_values in self.complex_values for name in complex_values}
        for name in names:
            column = np.full(n, float(specs[name].default))
            for k, complex_values in enumerate(self.complex_values):
                if name in complex_values:
                    column[k] = complex_values[name]
            values[name] = column
        bases, rows, ranks, nff, units = (np.array(column) for column in zip(*self.phytomers))
        for name, (positions, per_position) in ShootArchitecture._phytomer_values(dimensions, rows, ranks, nff,
                                                                                  units).items():
            column = values.setdefault(name, np.full(n, float(specs[name].default)))
            for position, array in zip(positions, per_position):
                column[bases + position] = array
        return NewNodes(parent=self.parent, edge=self.edge, label=self.element_label,
                        values=dict(element_index=self.element_index),
                        complexes=NewComplexes(anchor=self.anchor, scale=self.scale, label=self.label,
                                               values=values),
                        complex=self.complex_of)
