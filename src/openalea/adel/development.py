"""
A simple deterministic development for ShootArchitecture: ScenarioDevelopment writes what a development and a
senescence model write to it (the coupling contract, openalea.adel.architecture), from each plant's thermal time,
advanced at each call. It is a scenario, linear in thermal time, not Adel-Wheat's kinetics (which are not ported):

- per plant, its axis_number: its axes whose first leaf has appeared (the axes table's emf1);
- per axis, its phytomer_number: its metamers appeared, one per phyllochron from its first leaf, nf at most; its
  haun_stage: the phyllochrons since its first leaf, nf at most;
- per organ (durations in phyllochrons, from the metamer's appearance): the blade grows linearly to its final length
  (the phytomer table's Ll) in `blade_duration`; the sheath (Gl) from `sheath_start` for `sheath_duration`; the
  internode (El) from `internode_start` for `internode_duration`; the blade and the sheath are visible above the
  previous sheath, the internode for `internode_visible_fraction` of its length; the blade senesces from
  `senescence_start` for `senescence_duration`; its age is the thermal time since it appeared.

Its axes and organs are ShootArchitecture's, those with a row in the axes table (axis_row >= 0): on a plant another
model shares (a crown's axis and its subcrown internode, the roots' organs), the others keep their values.

A plant's thermal time starts at `start` (by default its main stem's first leaf) minus its `delay`. It is the
plant's at the start of each call: a call adds the degree-days of the step before it, the air temperature (the
forcing `air_temperature`, °C, read at that step's start) above `base_temperature` (none when that step started at
the call's time); without that forcing, each call adds `thermal_time_step`. A call grows the plants to their thermal
time whatever its step: a stand is grown to its stage by one call (its `start` the stage). ScenarioShoot is the plant
model of ShootArchitecture grown by it; shoot_scene builds a scene of such plants:

    scene = shoot_scene(stand(20, rows=4, plant_spacing=0.04, jitter=0.005, seed=0), delays=random_delays(20),
                        leaves=LeafShapes.echap(), thermal_time_step=10.)
    scene.run()                                        # 10 °Cd of development, then the shoot
    scene = shoot_scene(time_step=86400., forcings={"air_temperature": 15.}, start=0.)    # 15 °Cd a day
"""
from dataclasses import dataclass, replace

import numpy as np

from openalea.adel.architecture import SCALES, ShootArchitecture
from openalea.adel.tables import phytomer_dimensions, plant_rows, plant_tables, population_tables
from openalea.metafspm.coupling.component import (FunctionalComponent, forcing_variable, input_variable, parameter,
                                                   state_variable)
from openalea.metafspm.coupling.composite_wrapper import CompositeModel
from openalea.metafspm.scene.population import Planting, Stand, planting
from openalea.metafspm.scene.scene import Scene
from openalea.metafspm.solve.decorator import postsegmentation, potential

BY = "model_scenario_development"
ORGAN_KINDS = ("StemInternode", "Sheath", "Blade")
DAY = 86400.


@dataclass
class ScenarioDevelopment(FunctionalComponent):
    """The organs' lengths of a linear scenario, given at each call from the plants' thermal time."""

    # --- Read from ShootArchitecture ---
    shoot_tables: object = input_variable(unit='', description="The plants' tables (ShootArchitecture's)",
                                          by="ShootArchitecture", default=None, scale=SCALES.Plant, dtype='object')
    axis_row: int = input_variable(unit='adim', description="The axis' row in the population's axes table",
                                   by='model_shoot_architecture', default=-1, scale=SCALES.Axis, dtype='int')
    organ_rank: int = input_variable(unit='adim', description="The organ's phytomer rank",
                                     by='model_shoot_architecture', default=0, scale=SCALES.Organ, dtype='int')

    # --- Written to ShootArchitecture ---
    thermal_time: float = state_variable(unit='°Cd', description="The plant's thermal time in the scenario",
                                         value_comment="NaN: not started", kind="intensive",
                                         default=np.nan, scale=SCALES.Plant)
    axis_number: float = state_variable(unit='adim', description="The number of the plant's axes emerged",
                                        kind="descriptor", default=1., scale=SCALES.Plant)
    phytomer_number: float = state_variable(unit='adim', description="The number of phytomers the axis has",
                                            value_comment="1 at its creation: an axis appears with its first phytomer",
                                            kind="descriptor", default=1., scale=SCALES.Axis)
    haun_stage: float = state_variable(unit='adim',
                                       description="The axis' Haun stage: the phyllochrons since its first leaf",
                                       kind="descriptor", default=0., scale=SCALES.Axis)
    organ_length: float = state_variable(unit='m', description="The organ's length", kind="descriptor",
                                         default=0., scale=SCALES.Organ)
    visible_length: float = state_variable(unit='m', description="The organ's length out of the sheaths' tube",
                                           kind="descriptor", default=0., scale=SCALES.Organ)
    organ_senesced_length: float = state_variable(unit='m', description="The organ's senesced length",
                                                  kind="descriptor", default=0., scale=SCALES.Organ)
    age: float = state_variable(unit='°Cd', description="The blade's age", kind="descriptor",
                                default=0., scale=SCALES.Organ)

    # --- Forcing ---
    air_temperature: float = forcing_variable(unit='°C', description="The air temperature (its mean over the step)",
                                              required=False)

    # --- Parameters ---
    length_unit: float = parameter(unit='m', description="The tables' length unit in m", by=BY, default=0.01)
    phyllochron: float = parameter(unit='°Cd', description="The phyllochron of every axis", by=BY, default=88.92)
    thermal_time_step: float = parameter(unit='°Cd', description="The thermal time of a call",
                                         value_comment="used without the forcing air_temperature", by=BY, default=10.)
    base_temperature: float = parameter(unit='°C', description="The temperature above which the degree-days count",
                                        value_comment="with the forcing air_temperature", by=BY, default=0.)
    start: float = parameter(unit='°Cd', description="The scenario's start",
                             value_comment="NaN: the main stem's first leaf (emf1)", by=BY, default=np.nan)
    delay: float = parameter(unit='°Cd', description="The plant's delay: its scenario starts that much later", by=BY,
                             default=0.)
    blade_duration: float = parameter(unit='phyllochron',
                                      description="The blade's growth, from its phytomer's appearance", by=BY,
                                      default=1.6)
    sheath_start: float = parameter(unit='phyllochron',
                                    description="The sheath's growth start, after its phytomer's appearance", by=BY,
                                    default=1.)
    sheath_duration: float = parameter(unit='phyllochron', description="The sheath's growth", by=BY, default=1.)
    internode_start: float = parameter(unit='phyllochron',
                                       description="The internode's growth start, after its phytomer's appearance",
                                       by=BY, default=2.)
    internode_duration: float = parameter(unit='phyllochron', description="The internode's growth", by=BY, default=1.5)
    internode_visible_fraction: float = parameter(unit='adim', description="The internode's visible fraction", by=BY,
                                                  default=0.5)
    senescence_start: float = parameter(unit='phyllochron',
                                        description="The blade's senescence start, after its phytomer's appearance",
                                        by=BY, default=5.)
    senescence_duration: float = parameter(unit='phyllochron', description="The blade's senescence", by=BY, default=3.)

    @potential
    def _thermal_time(self, thermal_time, shoot_tables, start, delay, thermal_time_step, base_temperature):
        """The plants' thermal time advanced by the step before (set at start - delay at the first call); a plant's
        step, its parameters per plant."""
        thermal_time = np.asarray(thermal_time, dtype=float)
        first = self._start(self._tables(shoot_tables), start) - delay
        return np.where(np.isnan(thermal_time), first,
                        thermal_time + self._degree_days(thermal_time_step, base_temperature))

    @potential
    def _axis_number(self, thermal_time, shoot_tables):
        """Per plant, its axes emerged at its thermal time (in its tables' rows): before the shoot's segmentation,
        which creates the axes it asks for."""
        tables = self._tables(shoot_tables)
        thermal_time = np.asarray(thermal_time, dtype=float)
        return np.array([float(np.sum(thermal_time[p] >= tables["emf1"][tables["plant_rows"][number]]))
                         for p, number in enumerate(tables["numbers"].tolist())])

    @potential(location="Axis")
    def _axis_development(self, thermal_time, shoot_tables, axis_row, phytomer_number, haun_stage,
                          phyllochron) -> dict:
        """Per axis, at its plant's thermal time: its phytomers appeared and its Haun stage (an axis' step, its
        phyllochron its plant's); before the shoot's segmentation, which creates the phytomers it asks for."""
        return self._axes_given(thermal_time, shoot_tables, axis_row, phytomer_number, haun_stage, phyllochron)

    @potential(location="Organ")
    def _organ_development(self, thermal_time, shoot_tables, axis_row, organ_rank, organ_length, visible_length,
                           organ_senesced_length, age, phyllochron, length_unit, internode_visible_fraction,
                           blade_duration, sheath_start, sheath_duration, internode_start, internode_duration,
                           senescence_start, senescence_duration) -> dict:
        """Per organ, at its plant's thermal time: its lengths and a blade's age (an organ's step, its parameters
        its plant's)."""
        return self._organs_given(thermal_time, shoot_tables, axis_row, organ_rank, organ_length, visible_length,
                                  organ_senesced_length, age, phyllochron, length_unit, internode_visible_fraction,
                                  blade_duration, sheath_start, sheath_duration, internode_start, internode_duration,
                                  senescence_start, senescence_duration)

    def _degree_days(self, thermal_time_step, base_temperature):
        """Per plant, the thermal time a call adds: with the forcing air_temperature, the degree-days of the step
        before (its mean temperature over that step, above *base_temperature*, times its duration), kept from that
        step's call, none when that call started at the same time (no time passed: e.g. a stand grown to its stage
        by a call before a scene's first step); without it, *thermal_time_step*."""
        try:
            temperature = self.forcing("air_temperature")
        except KeyError:                                  # no air temperature given
            return np.asarray(thermal_time_step, dtype=float)
        step_before = self.kept("degree_days")
        start = self.forcing_time("start")
        above_base = np.maximum(0., temperature - np.asarray(base_temperature, dtype=float))
        self.keep("degree_days", (start, above_base * self.dt / DAY))
        return 0. if step_before is None or step_before[0] >= start else step_before[1]

    @postsegmentation(location="Axis")
    def _new_axes(self, thermal_time, shoot_tables, axis_row, phytomer_number, haun_stage, phyllochron) -> dict:
        """What is given per axis again, after the shoot's segmentation: the axes it created are given their
        phytomers and Haun stages at once (their tillers' inclination ramps with it), before the shoot's
        postsegmentation."""
        return self._axes_given(thermal_time, shoot_tables, axis_row, phytomer_number, haun_stage, phyllochron)

    @postsegmentation(location="Organ")
    def _new_organs(self, thermal_time, shoot_tables, axis_row, organ_rank, organ_length, visible_length,
                    organ_senesced_length, age, phyllochron, length_unit, internode_visible_fraction, blade_duration,
                    sheath_start, sheath_duration, internode_start, internode_duration, senescence_start,
                    senescence_duration) -> dict:
        """What is given per organ again, after the shoot's segmentation: the organs it created are given their
        lengths at once (as cnwgrass writes a new metamer's), before the shoot's postsegmentation."""
        return self._organs_given(thermal_time, shoot_tables, axis_row, organ_rank, organ_length, visible_length,
                                  organ_senesced_length, age, phyllochron, length_unit, internode_visible_fraction,
                                  blade_duration, sheath_start, sheath_duration, internode_start, internode_duration,
                                  senescence_start, senescence_duration)

    def _axes_given(self, thermal_time, shoot_tables, axis_row, phytomer_number, haun_stage, phyllochron):
        """{variable: values} per axis at its plant's *thermal_time*: its phytomer number (its phytomers appeared,
        one per *phyllochron* from its first leaf, nf at most) and its Haun stage; another model's axis (without a
        row, axis_row -1) keeps its *phytomer_number* and *haun_stage*."""
        tables = self._tables(shoot_tables)
        emf1, nf = tables["emf1"], tables["nf"]
        own = np.asarray(axis_row, dtype=int) >= 0
        axis_row = np.where(own, np.asarray(axis_row, dtype=int), 0)
        axis_time = np.asarray(thermal_time, dtype=float)[self.data_structure.owner("Plant", of="Axis")]
        phyllochron = np.asarray(phyllochron, dtype=float)
        ranks = np.arange(1, tables["Ll"].shape[1])
        appeared = emf1[axis_row][:, None] + (ranks[None, :] - 1) * phyllochron[:, None]
        appeared_number = np.sum((axis_time[:, None] >= appeared) & (ranks[None, :] <= nf[axis_row][:, None]),
                                 axis=1).astype(float)
        given_haun_stage = np.clip((axis_time - emf1[axis_row]) / phyllochron, 0., nf[axis_row])
        return dict(phytomer_number=np.where(own, appeared_number, np.asarray(phytomer_number, dtype=float)),
                    haun_stage=np.where(own, given_haun_stage, np.asarray(haun_stage, dtype=float)))

    def _organs_given(self, thermal_time, shoot_tables, axis_row, organ_rank, organ_length, visible_length,
                      organ_senesced_length, age, phyllochron, length_unit, internode_visible_fraction,
                      blade_duration, sheath_start, sheath_duration, internode_start, internode_duration,
                      senescence_start, senescence_duration):
        """{variable: values} per organ at its plant's *thermal_time*: its lengths (m), out of the sheaths' tube and
        senesced, and a blade's age; its parameters (durations in phyllochrons) one value per organ, its plant's.
        Another model's organ (on an axis without a row, or of another kind: a crown's, a root's) keeps its
        *organ_length*, *visible_length*, *organ_senesced_length* and *age*."""
        ds = self.data_structure
        tables = self._tables(shoot_tables)
        row = np.asarray(axis_row, dtype=int)[ds.owner("Axis", of="Organ")]
        kind = self._organ_kinds()
        own = (row >= 0) & (kind >= 0)
        row = np.where(own, row, 0)
        rank = np.asarray(organ_rank, dtype=int)
        time = np.asarray(thermal_time, dtype=float)[ds.owner("Plant", of="Organ")]
        phyllochron = np.asarray(phyllochron, dtype=float)
        timing = {name: np.asarray(values, dtype=float) * phyllochron for name, values in (
            ("blade_duration", blade_duration), ("sheath_start", sheath_start), ("sheath_duration", sheath_duration),
            ("internode_start", internode_start), ("internode_duration", internode_duration),
            ("senescence_start", senescence_start), ("senescence_duration", senescence_duration))}
        organ_unit = np.asarray(length_unit, dtype=float)
        organ_appeared = tables["emf1"][row] + (rank - 1) * phyllochron

        def ramp(begin, duration):
            return np.clip((time - begin) / duration, 0., 1.)

        def dimension(name, ranks_of):
            values = tables[name][row, np.clip(ranks_of, 0, None)]
            return np.where((ranks_of >= 1) & np.isfinite(values), values, 0.)

        blade = dimension("Ll", rank) * ramp(organ_appeared, timing["blade_duration"])
        sheath = dimension("Gl", rank) * ramp(organ_appeared + timing["sheath_start"], timing["sheath_duration"])
        internode = dimension("El", rank) * ramp(organ_appeared + timing["internode_start"],
                                                 timing["internode_duration"])
        previous_sheath = dimension("Gl", rank - 1) * ramp(organ_appeared - phyllochron + timing["sheath_start"],
                                                           timing["sheath_duration"])
        blade_senesced = np.minimum(dimension("Ll", rank) * ramp(organ_appeared + timing["senescence_start"],
                                                                 timing["senescence_duration"]) * (blade > 0), blade)
        visible_fraction = np.asarray(internode_visible_fraction, dtype=float)
        length = np.select([kind == 0, kind == 1, kind == 2], [internode, sheath, blade], 0.)
        visible = np.select([kind == 0, kind == 1, kind == 2],
                            [visible_fraction * internode, np.maximum(0., sheath - previous_sheath),
                             np.maximum(0., blade - previous_sheath)], 0.)
        senesced = np.where(kind == 2, blade_senesced, 0.)
        age = np.where(own & (kind == 2) & (time >= organ_appeared), time - organ_appeared,
                       np.asarray(age, dtype=float))              # a blade not appeared yet keeps its age
        return dict(organ_length=np.where(own, length * organ_unit, np.asarray(organ_length, dtype=float)),
                    visible_length=np.where(own, visible * organ_unit, np.asarray(visible_length, dtype=float)),
                    organ_senesced_length=np.where(own, senesced * organ_unit,
                                                   np.asarray(organ_senesced_length, dtype=float)), age=age)

    def _tables(self, shoot_tables):
        """From the population's tables (the plants' shoot_tables, ShootArchitecture's): emf1 and nf per row, each
        plant's number in them and their plants' rows, Ll, Gl, El per (row, rank) in their unit; read once per set of
        plants' tables."""
        records = tuple(shoot_tables)

        def build():
            axes, phytomers, numbers, _ = population_tables(records)
            dimensions = phytomer_dimensions(phytomers, axes, 1.)
            return dict(emf1=axes["emf1"].to_numpy(dtype=float), nf=axes["nf"].to_numpy(dtype=int), numbers=numbers,
                        plant_rows={int(number): plant_rows(axes, int(number)) for number in set(numbers.tolist())},
                        Ll=dimensions["Ll"], Gl=dimensions["Gl"], El=dimensions["El"])
        return self.cached("tables", build, key=tuple(id(record) for record in records))

    def _start(self, tables, start):
        """Each plant's start: its parameter *start*, by default (NaN) its main stem's first leaf (emf1)."""
        main_stem = np.array([tables["emf1"][tables["plant_rows"][number][0]] for number in tables["numbers"]])
        start = np.asarray(start, dtype=float)
        return np.where(np.isnan(start), main_stem, start)

    def _organ_kinds(self):
        """Per Organ: 0 internode, 1 sheath, 2 blade (from the MPG's labels), once per topology."""
        def build():
            kinds = {name: k for k, name in enumerate(ORGAN_KINDS)}
            return np.array([kinds.get(label, -1) for label in
                             self.data_structure.labels(location="Organ", per_node=False).tolist()], dtype=int)
        return self.cached("organ_kinds", build)


class ScenarioShoot(CompositeModel):
    """A plant model: ShootArchitecture grown by ScenarioDevelopment, a part other plant models nest (metafspm's
    CompositeModel parts: {"shoot": ScenarioShoot, ...}, its initiators by naming it in theirs). They need no
    translator: each input is the other's variable of the same name, shared on the DataStructure (the development's
    outputs are the shoot's inputs; the shoot's tables, axis rows and organ ranks the development's), as couplings()
    lists. Its one step, grow: both run together, in one schedule: the development's potential (its thermal time and
    what it gives per plant, axis and organ), the shoot's segmentation (the axes and phytomers asked for), the
    development's postsegmentation (the new axes' Haun stages, the new organs' lengths), then the shoot's
    postsegmentation (its elements). A call grows the plants to their thermal time, whatever its step."""
    initiators = (ShootArchitecture,)
    parts = {"development": ScenarioDevelopment, "architecture": ShootArchitecture}
    schedule = ("grow",)

    def grow(self):
        """The development and the shoot run together; then, while the axes the call created have phytomers to come
        (the development gives a new axis' phytomer_number in its postsegmentation, after the shoot's segmentation),
        the shoot's segmentation and the development's postsegmentation again, and the shoot's postsegmentation
        once more: a call creates every axis, phytomer and element its thermal time asks for, one large step as the
        small steps reaching it."""
        self.run_together(self.development, self.architecture)
        ds = self.architecture.data_structure
        grown = False
        while True:
            nodes = len(ds.entity_ids("node"))
            self.architecture.run_stage("segmentation")
            if len(ds.entity_ids("node")) == nodes:
                break
            self.development.run_stage("postsegmentation")
            grown = True
        if grown:
            self.architecture.run_stage("postsegmentation")


def random_delays(n_plants, max_delay=60., seed=0):
    """*n_plants* development delays (°Cd), uniform up to *max_delay*."""
    return np.random.default_rng(seed).uniform(0., max_delay, n_plants).tolist()


def shoot_scene(plants=None, model=ScenarioShoot, time_step="1 h", delays=None, plant_scenarios=None,
                depth="coleoptile", forcings=None, recorder=None, compute=None, **parameters):
    """A scene of *plants* (metafspm's Planting, m, or its stand(); by default one at the origin), each a *model*
    (ScenarioShoot or a model nesting it) with *parameters* (its components' parameters, e.g. axes_table,
    phytomer_table, leaves, n_sect, thermal_time_step, tiller_inclination, leaf_number), its own *plant_scenarios*
    (one dict per plant, e.g. its leaf_number) and *delays* (°Cd, ScenarioDevelopment's delay of each plant). A
    plant's shoot starts at its planting depth (metafspm's Planting.depth, below the soil surface, as a root model's
    would): *depth* (m) for every plant, "coleoptile" for each its coleoptile_sowing_depth, None for the Plantings'.
    *forcings*, *recorder*, *compute*: the Scene's (e.g. forcings={"air_temperature": ...}, °C, for
    ScenarioDevelopment's degree-days; compute its run configuration, metafspm's Compute)."""
    plants = plants or [Planting()]
    own = [dict(values) for values in (plant_scenarios or [{}] * len(plants))]
    for values, delay in zip(own, delays or [0.] * len(plants)):
        values["delay"] = float(delay)
    if depth is not None:
        depths = [coleoptile_sowing_depth(dict(parameters, **values)) if depth == "coleoptile" else float(depth)
                  for values in own]
        placed = [replace(plant, depth=value) for plant, value in zip(plants, depths)]
        plants = Stand(placed, plants.xrange, plants.yrange, plants.models) if isinstance(plants, Stand) else placed
    return Scene(planting(plants, model, parameters, own), time_step=time_step, forcings=forcings, recorder=recorder,
                 compute=compute)


COLEOPTILE_DEPTH_FRACTION = 0.9


def coleoptile_sowing_depth(scenario):
    """A plant's default sowing depth (m): COLEOPTILE_DEPTH_FRACTION of its coleoptile's final length (its main stem's
    rank-1 sheath, Gl, in its tables: tables.plant_tables), so that the coleoptile's tip emerges."""
    defaults = dict(table_plant=ShootArchitecture.table_plant, phyllochron=ShootArchitecture.phyllochron,
                    length_unit=ShootArchitecture.length_unit)
    axes, phytomers, number = plant_tables(scenario, defaults)
    first = phytomers.loc[(phytomers["plant"] == number) & (phytomers["axe"] == "MS") & (phytomers["n"] == 1), "Gl"]
    unit = float(scenario.get("length_unit", defaults["length_unit"]))
    return COLEOPTILE_DEPTH_FRACTION * float(first.iloc[0]) * unit if len(first) else 0.
