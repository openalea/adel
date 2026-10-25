"""
Adel's shoot architecture as a metafspm StructuralComponent: the plants' topology in the MPG (axes, phytomers, organs,
elements), the organs' dimensions and shapes, their elements, and their geometry for light models and display. The
organs' lengths are given by other models (development, senescence), as cnwgrass gives them to Adel's AdelDyn.

    table = pd.DataFrame([dict(plant="p0", model=MyModel, x=0., y=0., z=0., rotation=0.,
                               scenario={"parameters": {"axes_table": axeT, "phytomer_table": phytoT,
                                                        "table_plant": 1}})])

with MyModel.initiators = (ShootArchitecture,).

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

The tables are Adel's (its axeT and phytoT, one row per axis and per axis and rank, every plant of the population in
them, column 'plant'); `table_plant`, a per-plant parameter, says which of their plants a plant is.

The elements are placed by Adel's turtle (openalea.adel.turtle), from each plant's position and rotation: their base
x1, y1, z1, their end x2, y2, z2, and the frame they are drawn in (heading_x…, up_x…), at each step after their
lengths. triangles() gives their triangles (metafspm's coupling.geometry.Triangles) for light models and display.

The coupling contract, what other models write and read (lengths in m, on the DataStructure):

- a development model writes, per Plant, `axis_number` (its axes emerged, in its axes table's order); per Axis,
  `phytomer_number` (its phytomers initiated); per Organ, `organ_length`, `visible_length` (out of the sheaths' tube),
  a blade's `rolled_length` and `age`, and optionally a blade's `shape_mature_length` and `shape_max_width` (its
  final dimensions, the phytomer table's at its creation, as elongwheat's lamina_Lmax and leaf_Wmax);
- a senescence model writes, per Organ, `organ_senesced_length`;
- ShootArchitecture derives the elements from them at each call (their lengths, areas, green and senesced parts,
  placement): the elements are its own, the others read them, they do not write them;
- a growth or metabolism model reads the elements' `area`, `green_area`, `senesced_area`, the organs' `organ_area`;
- a light model reads triangles() (each triangle's element, its node index) and gives its results back per element
  (Triangles.per_owner).

The organs a call creates get their lengths at the next call; a development model giving them at once (as cnwgrass
after add_metamer) writes them after the call and calls ShootArchitecture again (comparison.given_development's
GivenShoot).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from openalea.metafspm.coupling.component import (NewComplexes, NewNodes, StructuralComponent, input_variable,
                                                   parameter, state_variable)
from openalea.metafspm.coupling.geometry import Triangles, place, slim_cylinders
from openalea.metafspm.data_structure.configs import PropsConfig, ScalesConfig
from openalea.metafspm.solve.decorator import postsegmentation, segmentation
from openalea.adel.blade_mesh import blade_element_triangles
from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.turtle import interpret, traversal_order

BY = "model_shoot_architecture"
DESCRIPTOR = "descriptor"
SCALES = ScalesConfig
ORGANS = ("internode", "sheath", "blade")


def doc(unit, description, value_comment=""):
    """The documentation entries of a declaration."""
    return dict(unit=unit, unit_comment="", description=description, min_value="", max_value="",
                value_comment=value_comment, references="", DOI="")


def outputs(k):
    """The return annotation of a step with k outputs besides its own: (values, name, values, ...)."""
    return tuple[(np.ndarray,) + (str, np.ndarray) * k]


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

    # --- Inputs: what development and senescence models write ---
    axis_number: float = input_variable(**doc('adim', "The number of the plant's axes emerged, in its axes table's "
                                                      "order (the main stem first)"),
                                        by='model_development', initialize=1., scale=SCALES.Plant)
    phytomer_number: float = input_variable(**doc('adim', "The number of phytomers the axis has (its development)"),
                                            by='model_development', initialize=0., scale=SCALES.Axis)
    organ_length: float = input_variable(**doc('m', "The organ's length"), by='model_development', initialize=0.,
                                         scale=SCALES.Organ)
    visible_length: float = input_variable(**doc('m', "The organ's length out of the sheaths' tube"),
                                           by='model_development', initialize=0., scale=SCALES.Organ)
    rolled_length: float = input_variable(**doc('m', "The blade's rolled length in the sheaths' tube"),
                                          by='model_development', initialize=0., scale=SCALES.Organ)
    organ_senesced_length: float = input_variable(**doc('m', "The organ's senesced length"),
                                                  by='model_senescence', initialize=0., scale=SCALES.Organ)
    age: float = input_variable(**doc('°Cd', "The blade's age (its shape class when the database has ages)"),
                                by='model_development', initialize=0., scale=SCALES.Organ)
    shape_mature_length: float = input_variable(**doc('m', "The blade's final length, its shape's length",
                                                      "at its creation, the phytomer table's Ll"),
                                                by='model_development', initialize=0., scale=SCALES.Organ)
    shape_max_width: float = input_variable(**doc('m', "The blade's maximal width, its shape's width",
                                                  "at its creation, the phytomer table's Lw"),
                                            by='model_development', initialize=0., scale=SCALES.Organ)

    # --- The axes (from the axes table) ---
    axis_row: int = state_variable(**doc('adim', "The axis' row in the axes table"), state_variable_type=DESCRIPTOR,
                                   initialize=-1, scale=SCALES.Axis, dtype='int')
    bearer_row: int = state_variable(**doc('adim', "The bearer axis' row in the axes table", "-1: the main stem"),
                                     state_variable_type=DESCRIPTOR, initialize=-1, scale=SCALES.Axis, dtype='int')
    insertion_rank: int = state_variable(**doc('adim', "The rank of the bearer's phytomer the axis is on"),
                                         state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.Axis, dtype='int')
    axis_phytomers: int = state_variable(**doc('adim', "The number of phytomers the axis has in the MPG"),
                                         state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.Axis, dtype='int')
    final_leaf_number: int = state_variable(**doc('adim', "The axis' final number of leaves (Adel's nf)"),
                                            state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.Axis,
                                            dtype='int')
    final_haun_stage: float = state_variable(**doc('adim', "The axis' final Haun stage (Adel's HS_final)"),
                                             state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.Axis)
    has_ear: int = state_variable(**doc('adim', "1 when the axis bears an ear"), state_variable_type=DESCRIPTOR,
                                  initialize=0, scale=SCALES.Axis, dtype='int')
    axis_azimuth: float = state_variable(**doc('°', "The axis' azimuth at its insertion (Adel's az_insertion), used "
                                                    "when its first stem element is inclined"),
                                         state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.Axis)

    # --- The phytomers and organs (from the phytomer table, at their creation) ---
    phytomer_rank: int = state_variable(**doc('adim', "The phytomer's rank on its axis, from 1"),
                                        state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.Phytomer,
                                        dtype='int')
    organ_rank: int = state_variable(**doc('adim', "The organ's phytomer rank"), state_variable_type=DESCRIPTOR,
                                     initialize=0, scale=SCALES.Organ, dtype='int')
    ntop: int = state_variable(**doc('adim', "The organ's rank from the axis' top leaf (1: the flag leaf)"),
                               state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.Organ, dtype='int')
    diameter: float = state_variable(**doc('m', "The stem organ's diameter (phytomer table: Ed, Gd)"),
                                     state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.Organ)
    organ_azimuth: float = state_variable(**doc('°', "The organ's azimuth (the internode's: phytomer table's Azim)"),
                                          state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.Organ)
    inclination: float = state_variable(**doc('°', "The organ's inclination"), state_variable_type=DESCRIPTOR,
                                        initialize=0., scale=SCALES.Organ)
    shape_class: int = state_variable(**doc('adim', "The blade's shape class in the leaf database (Lindex)",
                                            "-1: none"), state_variable_type=DESCRIPTOR, initialize=-1,
                                      scale=SCALES.Organ, dtype='int')
    shape_seed: int = state_variable(**doc('adim', "The blade's shape in its class (Lseed)", "-1: none"),
                                     state_variable_type=DESCRIPTOR, initialize=-1, scale=SCALES.Organ, dtype='int')

    organ_area: float = state_variable(**doc('m2', "The organ's visible area (its visible elements')"),
                                       state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.Organ)

    # --- The elements ---
    element_index: int = state_variable(**doc('adim', "The element's position in its organ, from 0"),
                                        state_variable_type=DESCRIPTOR, initialize=0, scale=SCALES.SubOrgan,
                                        dtype='int')
    length: float = state_variable(**doc('m', "The element's length"), state_variable_type=DESCRIPTOR,
                                   initialize=0., scale=SCALES.SubOrgan)
    area: float = state_variable(**doc('m2', "The element's area"), state_variable_type=DESCRIPTOR, initialize=0.,
                                 scale=SCALES.SubOrgan)
    green_length: float = state_variable(**doc('m', "The element's green length"), state_variable_type=DESCRIPTOR,
                                         initialize=0., scale=SCALES.SubOrgan)
    green_area: float = state_variable(**doc('m2', "The element's green area"), state_variable_type=DESCRIPTOR,
                                       initialize=0., scale=SCALES.SubOrgan)
    senesced_length: float = state_variable(**doc('m', "The element's senesced length"),
                                            state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    senesced_area: float = state_variable(**doc('m2', "The element's senesced area"), state_variable_type=DESCRIPTOR,
                                          initialize=0., scale=SCALES.SubOrgan)
    is_green: int = state_variable(**doc('adim', "1 when the element is more green than senesced"),
                                   state_variable_type=DESCRIPTOR, initialize=1, scale=SCALES.SubOrgan, dtype='int')
    rolled_part: float = state_variable(**doc('m', "The blade element's rolled length"),
                                        state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    rolled_diameter: float = state_variable(**doc('m', "The blade element's rolled part's diameter"),
                                            state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    width: float = state_variable(**doc('m', "The blade element's mean width"), state_variable_type=DESCRIPTOR,
                                  initialize=0., scale=SCALES.SubOrgan)
    flat_base: float = state_variable(**doc('adim', "The blade element's flat part's base, relative to the "
                                                    "blade's visible length (Adel's srb)", "NaN: no flat part"),
                                      state_variable_type=DESCRIPTOR, initialize=np.nan, scale=SCALES.SubOrgan)
    flat_top: float = state_variable(**doc('adim', "The blade element's flat part's top, relative to the blade's "
                                                   "visible length (Adel's srt)", "NaN: no flat part"),
                                     state_variable_type=DESCRIPTOR, initialize=np.nan, scale=SCALES.SubOrgan)

    # --- The elements' placement (the turtle's frame where each is drawn, and its end) ---
    x1: float = state_variable(**doc('m', "The element's base, x"), state_variable_type=DESCRIPTOR, initialize=0.,
                               scale=SCALES.SubOrgan)
    y1: float = state_variable(**doc('m', "The element's base, y"), state_variable_type=DESCRIPTOR, initialize=0.,
                               scale=SCALES.SubOrgan)
    z1: float = state_variable(**doc('m', "The element's base, z"), state_variable_type=DESCRIPTOR, initialize=0.,
                               scale=SCALES.SubOrgan)
    x2: float = state_variable(**doc('m', "The element's end (a blade element's: its rolled part's), x"),
                               state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    y2: float = state_variable(**doc('m', "The element's end (a blade element's: its rolled part's), y"),
                               state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    z2: float = state_variable(**doc('m', "The element's end (a blade element's: its rolled part's), z"),
                               state_variable_type=DESCRIPTOR, initialize=0., scale=SCALES.SubOrgan)
    heading_x: float = state_variable(**doc('adim', "The element's heading, x"), state_variable_type=DESCRIPTOR,
                                      initialize=0., scale=SCALES.SubOrgan)
    heading_y: float = state_variable(**doc('adim', "The element's heading, y"), state_variable_type=DESCRIPTOR,
                                      initialize=0., scale=SCALES.SubOrgan)
    heading_z: float = state_variable(**doc('adim', "The element's heading, z"), state_variable_type=DESCRIPTOR,
                                      initialize=1., scale=SCALES.SubOrgan)
    up_x: float = state_variable(**doc('adim', "The element's up, x (its blade's side: the shape's x is along it)"),
                                 state_variable_type=DESCRIPTOR, initialize=-1., scale=SCALES.SubOrgan)
    up_y: float = state_variable(**doc('adim', "The element's up, y"), state_variable_type=DESCRIPTOR, initialize=0.,
                                 scale=SCALES.SubOrgan)
    up_z: float = state_variable(**doc('adim', "The element's up, z"), state_variable_type=DESCRIPTOR, initialize=0.,
                                 scale=SCALES.SubOrgan)

    # --- Parameters ---
    axes_table: object = parameter(**doc('', "Adel's axeT: one row per axis of the population's plants (plant, axe, "
                                             "nf, HS_final, hasEar, ...)"), by=BY, default=None)
    phytomer_table: object = parameter(**doc('', "Adel's phytoT: one row per plant, axis and rank (plant, axe, n, Ll, "
                                                 "Lw, Gl, Gd, El, Ed, Azim, Lindex, Lseed)"), by=BY, default=None)
    table_plant: int = parameter(**doc('adim', "The plant's number in the tables (their column 'plant')"), by=BY,
                                 default=1)
    n_sect: int = parameter(**doc('adim', "The number of elements of a blade"), by=BY, default=1)
    length_unit: float = parameter(**doc('m', "The tables' length unit in m (Adel's cm: 0.01)"), by=BY, default=0.01)
    leaves: object = parameter(**doc('', "The leaf shape database (LeafShapes; None: Adel's samples)"), by=BY,
                               default=None)
    min_length: float = parameter(**doc('m', "The shortest element drawn (Adel's min_length, 0.01 cm); a blade's "
                                             "triangles smaller than its square / 100 are left out"), by=BY,
                                  default=1e-4)
    phyllochron: float = parameter(**doc('°Cd', "The main stem's phyllochron, for its blades' age classes in a "
                                                "dynamic leaf database"), by=BY, default=88.92)
    tiller_phyllochron: float = parameter(**doc('°Cd', "The tillers' phyllochron, for their blades' age classes in a "
                                                       "dynamic leaf database (Adel's T1's for every tiller)"), by=BY,
                                          default=91.6)

    # --- The plant's structure ---

    @classmethod
    def initiate_plant(cls, g, plant, parameters):
        """The plant's main stem (its Axis, GrowthUnit) and its first phytomer."""
        axes, phytomers = parameters["axes_table"], parameters["phytomer_table"]
        number = int(parameters.get("table_plant", cls.table_plant))
        n_sect = int(parameters.get("n_sect", cls.n_sect))
        unit = float(parameters.get("length_unit", cls.length_unit))
        rows = axes.reset_index(drop=True)
        mine = cls._plant_rows(rows, number)
        main_stem = mine[0]
        scales, labels = g.scales, g.labels
        axis = g.add_component(plant, **PropsConfig(scale=scales.Axis, edge_type='/', label=labels.Axis.Shoot))
        unit_vid = g.add_component(axis, **PropsConfig(scale=scales.GrowthUnit, edge_type='/',
                                                       label=labels.GrowthUnit.Shoot))
        for name, value in dict(cls._axis_values(rows, main_stem), phytomer_number=1., axis_phytomers=1).items():
            g.property(name)[axis] = value
        g.property("axis_number")[plant] = 1.
        dimensions = cls._dimensions(phytomers, rows, unit)
        complexes, elements = phytomer_layout(labels, n_sect)
        created = []
        for anchor, scale, label in complexes:
            parent = unit_vid if anchor == "unit" else created[anchor]
            created.append(g.add_component(parent, **PropsConfig(scale=scale, edge_type='/', label=label)))
        cls._set_phytomer(g, created, dimensions, main_stem, 1, int(rows.loc[main_stem, "nf"]))
        vids = []
        for complex_index, parent, edge, label, _, index in elements:
            if parent == "previous":
                vid = g.add_component(created[complex_index], **PropsConfig(scale=scales.SubOrgan, edge_type='/',
                                                                            label=label))
            else:
                vid = g.add_component_with_topo(created[complex_index], vids[parent], **PropsConfig(
                    scale=scales.SubOrgan, edge_type=edge, label=label))
            g.property("element_index")[vid] = index
            vids.append(vid)

    @staticmethod
    def _plant_rows(rows, number):
        """The rows of plant *number* in the axes table: its main stem first, then its tillers in table order."""
        mine = rows.index[rows["plant"] == number].tolist()
        main = [row for row in mine if rows.loc[row, "axe"] == "MS"]
        if not main:
            raise ValueError(f"no main stem of plant {number} in the axes table")
        return main + [row for row in mine if row not in main]

    @staticmethod
    def _axis_values(rows, row):
        """An axis' values from its row of the axes table (its bearer's row among the same plant's)."""
        name = rows.loc[row, "axe"]
        bearer, rank = bearer_of(name)
        same_plant = rows.index[(rows["plant"] == rows.loc[row, "plant"]) & (rows["axe"] == bearer)].tolist()
        haun_stage = rows.loc[row, "HS_final"]
        return dict(axis_row=row, bearer_row=same_plant[0] if same_plant else -1, insertion_rank=rank,
                    final_leaf_number=int(rows.loc[row, "nf"]), has_ear=int(rows.loc[row, "hasEar"]),
                    final_haun_stage=float(haun_stage) if pd.notna(haun_stage) else float(rows.loc[row, "nf"]))

    @staticmethod
    def _dimensions(phytomers, axes, unit):
        """{column: array (axis row, rank)} of the phytomer table, lengths in m; NaN where an axis has no rank."""
        columns = {"Ll": unit, "Lw": unit, "Gl": unit, "Gd": unit, "El": unit, "Ed": unit, "Azim": 1., "Lindex": 1.,
                   "Lseed": 1.}
        ranks = int(phytomers["n"].max())
        out = {name: np.full((len(axes), ranks + 1), np.nan) for name in columns}
        keyed = phytomers.set_index(["plant", "axe", "n"])
        for row, (plant, axe) in enumerate(zip(axes["plant"], axes["axe"])):
            if (plant, axe) not in {k[:2] for k in keyed.index}:
                continue
            these = keyed.loc[(plant, axe)]
            for name, factor in columns.items():
                out[name][row, these.index.to_numpy(dtype=int)] = these[name].to_numpy(dtype=float) * factor
        return out

    @classmethod
    def _set_phytomer(cls, g, created, dimensions, row, rank, nff):
        """A phytomer's and its organs' values on the MPG (initiate_plant's)."""
        for name, value in cls._phytomer_values(dimensions, np.array([row]), np.array([rank]),
                                                np.array([nff])).items():
            level, values = value
            for vid, v in zip([created[i] for i in level], values):
                g.property(name)[vid] = np.asarray(v).reshape(-1)[0].item()     # this one phytomer's

    @staticmethod
    def _phytomer_values(dimensions, rows, ranks, nff):
        """{variable: (complex positions in a phytomer's layout, values per phytomer per position)}: for phytomers
        of axis rows *rows* at ranks *ranks*; layout positions 0 phytomer, 1 internode, 2 sheath, 3 blade."""
        def at(name):
            return dimensions[name][rows, ranks]
        shape_class, shape_seed = at("Lindex"), at("Lseed")
        valid = (shape_class != -999) & (shape_seed != -999) & np.isfinite(shape_class) & np.isfinite(shape_seed)
        ntop = nff - ranks + 1
        return {
            "phytomer_rank": ([0], [ranks]),
            "organ_rank": ([1, 2, 3], [ranks, ranks, ranks]),
            "ntop": ([1, 2, 3], [ntop, ntop, ntop]),
            "diameter": ([1, 2], [at("Ed"), at("Gd")]),
            "organ_azimuth": ([1], [at("Azim")]),
            "shape_mature_length": ([3], [at("Ll")]),
            "shape_max_width": ([3], [at("Lw")]),
            "shape_class": ([3], [np.where(valid, shape_class, -1).astype(int)]),
            "shape_seed": ([3], [np.where(valid, shape_seed, -1).astype(int)]),
        }

    def _dimension_arrays(self):
        cache = self.__dict__.get("_dimensions_cache")
        if cache is None:
            axes = self.axes_table.reset_index(drop=True)
            cache = self._dimensions(self.phytomer_table, axes, float(np.asarray(self._parameter("length_unit"))[0]))
            self.__dict__["_dimensions_cache"] = cache
        return cache

    def _parameter(self, name):
        return self.parameter_values(name)

    def _axes_lookup(self):
        """(rows of each table plant, values of each row) of the axes table (_plant_rows, _axis_values), read once."""
        cache = self.__dict__.get("_axes_cache")
        if cache is None:
            rows = self.axes_table.reset_index(drop=True)
            plant_rows = {int(number): self._plant_rows(rows, int(number)) for number in set(rows["plant"])}
            cache = (plant_rows, {row: self._axis_values(rows, row) for row in range(len(rows))})
            self.__dict__["_axes_cache"] = cache
        return cache

    @segmentation
    def _axis_phytomers(self, axis_phytomers, phytomer_number, axis_row) -> tuple[np.ndarray, str, NewNodes]:
        """The axes the plants' axis_number asks for (in their axes table's order, each with its first phytomer, once
        its bearer has the phytomer it is inserted on), and the phytomers the axes' phytomer_number asks for: created
        with their organs and elements (NewNodes), their dimensions and shapes from the phytomer table."""
        ds = self.data_structure
        g = ds._mtg
        plant_rows, axis_values = self._axes_lookup()
        count = np.array(axis_phytomers, dtype=int)
        target = np.floor(np.asarray(phytomer_number, dtype=float)).astype(int)
        axes = np.asarray(ds.entity_ids("Axis"))
        axis_of = {(g.complex(int(axis)), int(row)): int(axis) for axis, row in zip(axes, axis_row)}
        batch = _Batch(g.labels, int(np.asarray(self.parameter_values("n_sect"))[0]))

        # New axes, with their first phytomer
        plants = np.asarray(ds.entity_ids("Plant"))
        wanted = np.floor(np.asarray(ds.get("axis_number"), dtype=float)).astype(int)
        table_plant = np.asarray(ds.parameter_view("table_plant", "Plant")).astype(int)
        for plant, number, table in zip(plants.tolist(), wanted.tolist(), table_plant.tolist()):
            mine = plant_rows[table]
            existing = [row for row in mine if (plant, row) in axis_of]
            for row in mine[len(existing):number]:
                values = axis_values[row]
                predecessor = self._predecessor(g, None, 1, axis_of, plant, values["bearer_row"],
                                                values["insertion_rank"])
                if predecessor is None:
                    break                                   # its bearer's phytomer is not there yet
                batch.axis(plant, dict(values, axis_phytomers=1))
                batch.phytomer(None, predecessor, row, 1, values["final_leaf_number"], branched=True)

        # New phytomers on the axes
        for a, axis in enumerate(axes.tolist()):
            values = axis_values[int(axis_row[a])]
            for rank in range(count[a] + 1, target[a] + 1):
                predecessor = self._predecessor(g, axis, rank, axis_of, g.complex(axis), values["bearer_row"],
                                                values["insertion_rank"])
                batch.phytomer(g.components(axis)[0], predecessor, int(axis_row[a]), rank,
                               values["final_leaf_number"], branched=False)
                count[a] = rank
        return count, "new_phytomers", batch.new_nodes(self._dimension_arrays(), self._variable_specs)

    # --- The elements ---

    def _kinds(self):
        """Per node: its element kind (0 hidden, 1 stem, 2 leaf) and its organ's (0 internode, 1 sheath, 2 blade),
        read from the MPG's labels once per topology."""
        ds = self.data_structure
        key = (ds.topology_version, ds.n_nodes())
        cache = self.__dict__.get("_kinds_cache")
        if cache is None or cache[0] != key:
            g = ds._mtg
            labels = g.labels
            element = {labels.SubOrgan.HiddenElement: 0, labels.SubOrgan.StemElement: 1,
                       labels.SubOrgan.LeafElement: 2}
            organ = {labels.Organ.StemInternode: 0, labels.Organ.Sheath: 1, labels.Organ.Blade: 2}
            nodes = ds.entity_ids("node").tolist()
            cache = (key, np.array([element.get(g.label(v), -1) for v in nodes]),
                     np.array([organ.get(g.label(g.complex(v)), -1) for v in nodes]))
            self.__dict__["_kinds_cache"] = cache
        return cache[1], cache[2]

    def _leaf_shapes(self):
        if self.leaves is None:
            self.__dict__.setdefault("_default_leaves", LeafShapes.samples())
            return self.__dict__["_default_leaves"]
        return self.leaves

    @postsegmentation
    def _length(self, length, element_index, organ_length, visible_length, organ_senesced_length, rolled_length,
                diameter, shape_mature_length, shape_max_width, shape_class) -> outputs(12):
        """The elements from their organs' lengths (Adel's internode_elements, sheath_elements, blade_elements):
        a stem organ's hidden element (its length out of sight) and visible one, its green and senesced parts; a
        blade's hidden part, and its sectors along its mature shape (its visible part the shape's top), each one's
        green and senesced parts, rolled part, flat part and width; areas from the diameters and the shapes' width
        profiles; the organs' visible areas."""
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
        sector = blade & (element_kind == 2) & has_shape & (l_visible > 1e-6)
        if sector.any():
            self._blade_sectors(out, is_green, sector, element_index, l_visible, l_senesced, l_rolled, mature,
                                max_width, classes, organ, shapes)

        # The organs' visible areas
        visible = ((element_kind == 1) & stem) | ((element_kind == 2) & blade)
        areas = np.bincount(organ[visible], weights=out["area"][visible], minlength=len(ds.entity_ids("Organ")))
        return (out["length"], "area", out["area"], "green_length", out["green_length"], "green_area",
                out["green_area"], "senesced_length", out["senesced_length"], "senesced_area", out["senesced_area"],
                "is_green", is_green, "rolled_part", out["rolled_part"], "rolled_diameter", out["rolled_diameter"],
                "width", out["width"], "flat_base", out["flat_base"], "flat_top", out["flat_top"], "organ_area",
                areas)

    @staticmethod
    def _blade_sectors(out, is_green, sector, element_index, l_visible, l_senesced, l_rolled, mature, max_width,
                       classes, organ, shapes):
        """Each visible blade's sectors (Adel's blade_elements): sector k covers [(k - 1), k] · L / n_sect of the
        mature shape (L its length), its visible part above L - visible length, its senesced part above L -
        senesced length, its flat part above L - (visible - rolled) length."""
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
        c = classes[sector]
        green_area = np.where(green_part > 0, shapes.relative_area(c, base / L, green_top / L) * L * w, 0.)
        senesced_area = np.where(senesced_part > 0, shapes.relative_area(c, green_top / L, top / L) * L * w, 0.)
        total = green_area + senesced_area
        flat_part = np.minimum(visible_part, np.maximum(0., top - limit_rolled))
        has_flat = flat_part > 0
        flat_bottom = top - flat_part
        flat_green_top = np.minimum(top, np.maximum(flat_bottom, limit_senesced))
        with np.errstate(divide="ignore", invalid="ignore"):
            rolled_part = np.where(flat_part < visible_part, visible_part - flat_part, 0.)
            rolled_diameter = np.where(rolled_part > 0, total * rolled_part / visible_part / np.pi / rolled_part, 0.)
            green_width = np.where(green_part > 0, shapes.mean_width(c, base / L, green_top / L) * w, 0.)
            senesced_width = np.where(senesced_part > 0, shapes.mean_width(c, green_top / L, top / L) * w, 0.)
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
        """Per node: its position in Adel's traversal (turtle.traversal_order), and whether it is its plant's base;
        once per topology."""
        ds = self.data_structure
        key = (ds.topology_version, ds.n_nodes())
        cache = self.__dict__.get("_traversal_cache")
        if cache is None or cache[0] != key:
            g = ds._mtg
            nodes = ds.entity_ids("node").tolist()
            base = np.array([g.parent(v) is None for v in nodes])
            cache = (key, traversal_order(g, nodes), base)
            self.__dict__["_traversal_cache"] = cache
        return cache[1], cache[2]

    @postsegmentation
    def _x1(self, length, rolled_part, inclination, organ_azimuth, axis_azimuth, bearer_row,
            phytomer_rank) -> outputs(11):
        """The elements' placement, their frames by Adel's turtle (turtle.interpret): each one's base, heading and up
        where it is drawn, and its end; from its plant's position and rotation, on the elements' lengths."""
        ds = self.data_structure
        g = ds._mtg
        element_kind, organ_kind = self._kinds()
        order, plant_base = self._traversal()
        organ, axis = np.asarray(ds.owner("Organ")), np.asarray(ds.owner("Axis"))
        plant, phytomer = np.asarray(ds.owner("Plant")), np.asarray(ds.owner("Phytomer"))
        plants = ds.entity_ids("Plant").tolist()
        plant_position = np.array([[g.property(name).get(v, 0.) for name in ("x", "y", "z")] for v in plants],
                                  dtype=float).reshape(-1, 3)
        plant_rotation = np.array([g.property("rotation").get(v, 0.) for v in plants], dtype=float)
        axis_is_main = np.asarray(bearer_row) < 0
        stem = element_kind == 1
        position, heading, up, end = interpret(
            order, element_kind.astype(np.int64), plant_base, axis.astype(np.int64), axis_is_main,
            (np.asarray(phytomer_rank)[phytomer] == 1) & ~axis_is_main[axis], plant_position[plant],
            plant_rotation[plant], np.asarray(axis_azimuth, dtype=float)[axis],
            np.where(stem, np.asarray(inclination, dtype=float)[organ], 0.),
            np.where(stem & (organ_kind == 0), np.asarray(organ_azimuth, dtype=float)[organ], 0.),
            np.asarray(length, dtype=float), np.asarray(rolled_part, dtype=float))
        return (position[:, 0], "y1", position[:, 1], "z1", position[:, 2], "x2", end[:, 0], "y2", end[:, 1], "z2",
                end[:, 2], "heading_x", heading[:, 0], "heading_y", heading[:, 1], "heading_z", heading[:, 2],
                "up_x", up[:, 0], "up_y", up[:, 1], "up_z", up[:, 2])

    def frames(self):
        """(position, heading, up), each (n, 3): the elements' frames where they are drawn."""
        ds = self.data_structure
        names = (("x1", "y1", "z1"), ("heading_x", "heading_y", "heading_z"), ("up_x", "up_y", "up_z"))
        return tuple(np.stack([np.asarray(ds.get(name), dtype=float) for name in frame], axis=1) for frame in names)

    def triangles(self, where=None) -> Triangles:
        """The elements' triangles (metafspm's coupling.geometry.Triangles), as Adel draws them: a stem element's
        slim cylinder of its organ's diameter; a blade element's rolled part's, of its rolled diameter, then its
        flat part's mesh above it (blade_mesh). *where*: a node mask."""
        ds = self.data_structure
        element_kind, _ = self._kinds()
        organ = np.asarray(ds.owner("Organ"))
        length, rolled = np.asarray(ds.get("length")), np.asarray(ds.get("rolled_part"))
        minimum = float(np.asarray(self.parameter_values("min_length"))[0])
        selected = np.ones(length.size, dtype=bool) if where is None else np.asarray(where, dtype=bool)
        leaf = selected & (element_kind == 2) & (length > 0)
        stem = selected & (element_kind == 1) & (length > 0) & (length >= minimum)
        cylinder = np.flatnonzero(stem | (leaf & (rolled > 0)))
        cylinder_length = np.where(stem, length, rolled)[cylinder]
        radius = np.where(stem, 0.5 * np.asarray(ds.get("diameter"))[organ],
                          0.5 * np.asarray(ds.get("rolled_diameter")))[cylinder]
        position, heading, up = self.frames()
        vertices = place(slim_cylinders(cylinder_length, radius, radius), position[cylinder], heading[cylinder],
                         up[cylinder]).reshape(-1, 3, 3)
        owner = np.repeat(cylinder, 8)
        blade_vertices, blade_owner = self._blade_triangles(leaf, rolled, minimum, organ)
        if blade_owner.size:
            placed = place(blade_vertices[:, None], position[blade_owner], heading[blade_owner],
                           up[blade_owner])[:, 0]
            vertices, owner = np.concatenate([vertices, placed]), np.concatenate([owner, blade_owner])
        return Triangles(vertices, owner)

    def _blade_triangles(self, leaf, rolled, minimum, organ):
        """The flat parts' meshes of the blade elements *leaf* (a node mask) in their local frames, above their rolled
        parts: (T, 3, 3) vertices, (T,) node indices. Each blade's shape is its class', seed's and, in a dynamic
        database, its age class' (its age in its axis' phyllochrons, minus 0.3, as Adel's). An element's mesh is
        kept from one call to the next while what it is made from is the same (its shape, visible length, flat
        part, ...): only the grown and new ones are meshed again, every one placed again."""
        ds = self.data_structure
        shapes = self._leaf_shapes()
        visible = np.asarray(ds.get("visible_length"))[organ]
        flat_base, flat_top = np.asarray(ds.get("flat_base")), np.asarray(ds.get("flat_top"))
        meshed = np.flatnonzero(leaf & (visible >= minimum) & np.isfinite(flat_base))
        if meshed.size == 0 or not shapes.fitted:
            return np.zeros((0, 3, 3)), np.zeros(0, dtype=np.int64)
        blade = organ[meshed]
        main = ds.map_to(np.asarray(ds.get("bearer_row"), dtype=float), "Axis", "Organ", "broadcast")[blade] < 0
        phyllochron = np.where(main, float(np.asarray(self.parameter_values("phyllochron"))[0]),
                               float(np.asarray(self.parameter_values("tiller_phyllochron"))[0]))
        age_class = shapes.age_class(np.asarray(ds.get("age"))[blade] / phyllochron - 0.3)
        shape = shapes.shape_index(np.asarray(ds.get("shape_class"))[blade], np.asarray(ds.get("shape_seed"))[blade],
                                   age_class)
        inclination = 1. if shapes.dynamic else None          # a dynamic database's shapes carry their inclination
        mature, max_width = np.asarray(ds.get("shape_mature_length")), np.asarray(ds.get("shape_max_width"))
        blade_inclination = np.asarray(ds.get("inclination"))
        min_area = minimum ** 2 / 100.
        vertex_ids = ds.entity_ids("node")
        previous = self.__dict__.get("_blade_mesh_cache", {})
        cache = {}
        parts, owners = [], []
        for node, organ_index, shape_index in zip(meshed.tolist(), blade.tolist(), shape.tolist()):
            if shape_index < 0:
                continue
            made_from = (id(shapes), shape_index,
                         blade_inclination[organ_index] if inclination is None else inclination, mature[organ_index],
                         max_width[organ_index], visible[node], flat_base[node], flat_top[node], rolled[node],
                         min_area)
            vertex = int(vertex_ids[node])
            kept = previous.get(vertex)
            if kept is not None and kept[0] == made_from:
                triangles = kept[1]
            else:
                triangles = blade_element_triangles(*shapes.fitted[shape_index], *made_from[2:8], min_area)
                triangles[..., 2] += rolled[node]
            cache[vertex] = (made_from, triangles)
            parts.append(triangles)
            owners.append(np.full(len(triangles), node, dtype=np.int64))
        self.__dict__["_blade_mesh_cache"] = cache
        if not parts:
            return np.zeros((0, 3, 3)), np.zeros(0, dtype=np.int64)
        return np.concatenate(parts), np.concatenate(owners)

    @staticmethod
    def _predecessor(g, axis, rank, axis_of, plant, bearer_row, insertion_rank):
        """The element a new phytomer's internode follows: the previous internode's last element; for an axis' first
        phytomer, its bearer's phytomer's internode's first element (NO_PARENT for a main stem); None while that
        phytomer does not exist."""
        if rank > 1:
            phytomer = g.components(g.components(axis)[0])[rank - 2]
        elif bearer_row < 0:
            return NewNodes.NO_PARENT
        else:
            bearer = axis_of.get((plant, bearer_row))
            if bearer is None:
                return None
            phytomers = g.components(g.components(bearer)[0])
            if len(phytomers) < max(insertion_rank, 1):
                return None
            return g.components(g.components(phytomers[max(insertion_rank, 1) - 1])[0])[0]
        internode = g.components(phytomer)[0]
        return g.components(internode)[-1]


class _Batch:
    """The complexes and elements of a step's new axes and phytomers, gathered for one NewNodes."""

    def __init__(self, labels, n_sect):
        self.labels = labels
        self.complexes_layout, self.elements_layout = phytomer_layout(labels, n_sect)
        self.anchor, self.scale, self.label, self.complex_values = [], [], [], []
        self.parent, self.edge, self.element_label, self.complex_of, self.element_index = [], [], [], [], []
        self.phytomers = []                                 # (layout's first complex, axis row, rank, nff)
        self.last_unit = None

    def _complex(self, anchor, scale, label, values=None):
        self.anchor.append(anchor)
        self.scale.append(scale)
        self.label.append(label)
        self.complex_values.append(values or {})
        return len(self.anchor) - 1

    def axis(self, plant, values):
        """A new Axis of *plant* with its GrowthUnit; the next phytomer is its first."""
        axis = self._complex(plant, SCALES.Axis, self.labels.Axis.Shoot, values)
        self.last_unit = self._complex(-1 - axis, SCALES.GrowthUnit, self.labels.GrowthUnit.Shoot)

    def phytomer(self, unit, predecessor, row, rank, nff, branched):
        """A phytomer under the GrowthUnit *unit* (None: the new axis' just added), its internode after
        *predecessor* (an element's id, or NO_PARENT)."""
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
        self.phytomers.append((base, row, rank, nff))

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
        bases, rows, ranks, nff = (np.array(column) for column in zip(*self.phytomers))
        for name, (positions, per_position) in ShootArchitecture._phytomer_values(dimensions, rows, ranks,
                                                                                  nff).items():
            column = values.setdefault(name, np.full(n, float(specs[name].default)))
            for position, array in zip(positions, per_position):
                column[bases + position] = array
        return NewNodes(parent=self.parent, edge=self.edge, label=self.element_label,
                        values=dict(element_index=self.element_index),
                        complexes=NewComplexes(anchor=self.anchor, scale=self.scale, label=self.label,
                                               values=values),
                        complex=self.complex_of)
