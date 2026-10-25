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

A plant's thermal time starts at `start` (by default its main stem's first leaf) minus its `delay`. ScenarioShoot is
the plant model of ShootArchitecture grown by it; shoot_scene builds a scene of such plants:

    scene = shoot_scene(sown_stand(20), leaves=LeafShapes.echap(), thermal_time_step=10.)
    scene.run()                                        # 10 °Cd of development, then the shoot
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from openalea.adel.architecture import SCALES, ShootArchitecture, doc
from openalea.adel.tables import population_tables
from openalea.metafspm.coupling.component import FunctionalComponent, input_variable, parameter, state_variable
from openalea.metafspm.scene.scene import Scene
from openalea.metafspm.solve.decorator import state

BY = "model_scenario_development"
ORGAN_KINDS = ("StemInternode", "Sheath", "Blade")


@dataclass
class ScenarioDevelopment(FunctionalComponent):
    """The organs' lengths of a linear scenario, given at each call from the plants' thermal time."""

    # --- Read from ShootArchitecture ---
    axis_row: int = input_variable(**doc('adim', "The axis' row in the population's axes table"),
                                   by='model_shoot_architecture', initialize=-1, scale=SCALES.Axis, dtype='int')
    organ_rank: int = input_variable(**doc('adim', "The organ's phytomer rank"), by='model_shoot_architecture',
                                     initialize=0, scale=SCALES.Organ, dtype='int')

    # --- Written to ShootArchitecture ---
    thermal_time: float = state_variable(**doc('°Cd', "The plant's thermal time in the scenario", "NaN: not started"),
                                         state_variable_type="intensive", initialize=np.nan, scale=SCALES.Plant)
    axis_number: float = state_variable(**doc('adim', "The number of the plant's axes emerged"),
                                        state_variable_type="descriptor", initialize=1., scale=SCALES.Plant)
    phytomer_number: float = state_variable(**doc('adim', "The number of phytomers the axis has"),
                                            state_variable_type="descriptor", initialize=0., scale=SCALES.Axis)
    haun_stage: float = state_variable(**doc('adim', "The axis' Haun stage: the phyllochrons since its first leaf"),
                                       state_variable_type="descriptor", initialize=0., scale=SCALES.Axis)
    organ_length: float = state_variable(**doc('m', "The organ's length"), state_variable_type="descriptor",
                                         initialize=0., scale=SCALES.Organ)
    visible_length: float = state_variable(**doc('m', "The organ's length out of the sheaths' tube"),
                                           state_variable_type="descriptor", initialize=0., scale=SCALES.Organ)
    organ_senesced_length: float = state_variable(**doc('m', "The organ's senesced length"),
                                                  state_variable_type="descriptor", initialize=0., scale=SCALES.Organ)
    age: float = state_variable(**doc('°Cd', "The blade's age"), state_variable_type="descriptor", initialize=0.,
                                scale=SCALES.Organ)

    # --- Parameters ---
    length_unit: float = parameter(**doc('m', "The tables' length unit in m"), by=BY, default=0.01)
    phyllochron: float = parameter(**doc('°Cd', "The phyllochron of every axis"), by=BY, default=88.92)
    thermal_time_step: float = parameter(**doc('°Cd', "The thermal time of a call"), by=BY, default=10.)
    start: float = parameter(**doc('°Cd', "The scenario's start", "NaN: the main stem's first leaf (emf1)"), by=BY,
                             default=np.nan)
    delay: float = parameter(**doc('°Cd', "The plant's delay: its scenario starts that much later"), by=BY,
                             default=0.)
    blade_duration: float = parameter(**doc('phyllochron', "The blade's growth, from its metamer's appearance"),
                                      by=BY, default=1.6)
    sheath_start: float = parameter(**doc('phyllochron', "The sheath's growth start, after its metamer's appearance"),
                                    by=BY, default=1.)
    sheath_duration: float = parameter(**doc('phyllochron', "The sheath's growth"), by=BY, default=1.)
    internode_start: float = parameter(**doc('phyllochron', "The internode's growth start, after its metamer's "
                                                            "appearance"), by=BY, default=2.)
    internode_duration: float = parameter(**doc('phyllochron', "The internode's growth"), by=BY, default=1.5)
    internode_visible_fraction: float = parameter(**doc('adim', "The internode's visible fraction"), by=BY,
                                                  default=0.5)
    senescence_start: float = parameter(**doc('phyllochron', "The blade's senescence start, after its metamer's "
                                                             "appearance"), by=BY, default=5.)
    senescence_duration: float = parameter(**doc('phyllochron', "The blade's senescence"), by=BY, default=3.)

    @state
    def _thermal_time(self, thermal_time) -> tuple[np.ndarray, str, np.ndarray, str, np.ndarray, str, np.ndarray,
                                                   str, np.ndarray, str, np.ndarray, str, np.ndarray, str,
                                                   np.ndarray]:
        """The plants' thermal time advanced by a step (set at start - delay at the first call), and what is given
        at it."""
        thermal_time = np.asarray(thermal_time, dtype=float)
        first = self._start() - self._plant_values("delay")
        thermal_time = np.where(np.isnan(thermal_time), first, thermal_time + self._plant_values("thermal_time_step"))
        given = self.given(thermal_time)
        return (thermal_time, "axis_number", given["axis_number"], "phytomer_number", given["phytomer_number"],
                "haun_stage", given["haun_stage"], "organ_length", given["organ_length"], "visible_length",
                given["visible_length"], "organ_senesced_length", given["organ_senesced_length"], "age",
                given["age"])

    def write(self):
        """Write what is given at the plants' current thermal time again: the organs created since are given their
        lengths (as cnwgrass writes a new metamer's at once)."""
        ds = self.data_structure
        for name, values in self.given(np.asarray(ds.get("thermal_time"), dtype=float)).items():
            ds.set(name, values)

    def given(self, thermal_time):
        """{variable: values} at the plants' *thermal_time*: their axis numbers, their axes' phytomer numbers and Haun
        stages, their organs' lengths (m) and blades' ages."""
        ds = self.data_structure
        tables = self._tables()
        emf1, nf = tables["emf1"], tables["nf"]

        # Plants: their axes emerged (in their tables' rows)
        axis_number = np.array([float(np.sum(thermal_time[p] >= emf1[tables["plant_rows"][number]]))
                                for p, number in enumerate(tables["numbers"].tolist())])

        # Axes: their metamers appeared
        axis_row = np.asarray(ds.get("axis_row"), dtype=int)
        axis_time = ds.map_to(thermal_time, "Plant", "Axis", "broadcast")
        axis_phyllochron = ds.map_to(self._plant_values("phyllochron"), "Plant", "Axis", "broadcast")
        ranks = np.arange(1, tables["Ll"].shape[1])
        appeared = emf1[axis_row][:, None] + (ranks[None, :] - 1) * axis_phyllochron[:, None]
        phytomer_number = np.sum((axis_time[:, None] >= appeared) & (ranks[None, :] <= nf[axis_row][:, None]),
                                 axis=1).astype(float)
        haun_stage = np.clip((axis_time - emf1[axis_row]) / axis_phyllochron, 0., nf[axis_row])

        # Organs: their lengths
        n_axes = len(ds.entity_ids("Axis"))
        row = axis_row[ds.map_to(np.arange(n_axes, dtype=float), "Axis", "Organ", "broadcast").astype(int)]
        rank = np.asarray(ds.get("organ_rank"), dtype=int)

        def per_organ(values):
            return ds.map_to(np.asarray(values, dtype=float), "Plant", "Organ", "broadcast")

        time = per_organ(thermal_time)
        phyllochron = per_organ(self._plant_values("phyllochron"))
        timing = {name: per_organ(self._plant_values(name)) * phyllochron for name in (
            "blade_duration", "sheath_start", "sheath_duration", "internode_start", "internode_duration",
            "senescence_start", "senescence_duration")}
        organ_unit = per_organ(self._plant_values("length_unit"))
        kind = self._organ_kinds()
        organ_appeared = emf1[row] + (rank - 1) * phyllochron

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
        visible_fraction = per_organ(self._plant_values("internode_visible_fraction"))
        length = np.select([kind == 0, kind == 1, kind == 2], [internode, sheath, blade], 0.)
        visible = np.select([kind == 0, kind == 1, kind == 2],
                            [visible_fraction * internode, np.maximum(0., sheath - previous_sheath),
                             np.maximum(0., blade - previous_sheath)], 0.)
        senesced = np.where(kind == 2, blade_senesced, 0.)
        age = np.where((kind == 2) & (time >= organ_appeared), time - organ_appeared,
                       np.asarray(ds.get("age"), dtype=float))              # a blade not appeared yet keeps its age
        return dict(axis_number=axis_number, phytomer_number=phytomer_number, haun_stage=haun_stage,
                    organ_length=length * organ_unit, visible_length=visible * organ_unit,
                    organ_senesced_length=senesced * organ_unit, age=age)

    def _plant_values(self, name):
        return np.asarray(self.data_structure.parameter_view(name, "Plant"), dtype=float)

    def _tables(self):
        """From the population's tables (ShootArchitecture's): emf1 and nf per row, each plant's number in them and
        their plants' rows, Ll, Gl, El per (row, rank) in their unit; read once per plant set."""
        ds = self.data_structure
        plants = tuple(ds.entity_ids("Plant").tolist())
        cache = self.__dict__.get("_tables_cache")
        if cache is None or cache[0] != plants:
            axes, phytomers, numbers, _ = population_tables(ds._mtg, plants)
            dimensions = ShootArchitecture._dimensions(phytomers, axes, 1.)
            tables = dict(emf1=axes["emf1"].to_numpy(dtype=float), nf=axes["nf"].to_numpy(dtype=int), numbers=numbers,
                          plant_rows={int(number): ShootArchitecture._plant_rows(axes, int(number))
                                      for number in set(numbers.tolist())},
                          Ll=dimensions["Ll"], Gl=dimensions["Gl"], El=dimensions["El"])
            cache = (plants, tables)
            self.__dict__["_tables_cache"] = cache
        return cache[1]

    def _start(self):
        """Each plant's start: its parameter, by default its main stem's first leaf (emf1)."""
        tables = self._tables()
        main_stem = np.array([tables["emf1"][tables["plant_rows"][number][0]] for number in tables["numbers"]])
        start = self._plant_values("start")
        return np.where(np.isnan(start), main_stem, start)

    def _organ_kinds(self):
        """Per Organ: 0 internode, 1 sheath, 2 blade (from the MPG's labels), once per topology."""
        ds = self.data_structure
        key = (ds.topology_version, len(ds.entity_ids("Organ")))
        cache = self.__dict__.get("_kinds_cache")
        if cache is None or cache[0] != key:
            g = ds._mtg
            kinds = {getattr(g.labels.Organ, name): k for k, name in enumerate(ORGAN_KINDS)}
            cache = (key, np.array([kinds.get(g.label(int(v)), -1) for v in ds.entity_ids("Organ")], dtype=int))
            self.__dict__["_kinds_cache"] = cache
        return cache[1]


class ScenarioShoot:
    """A plant model: ShootArchitecture grown by ScenarioDevelopment. At each call the development first, then the
    shoot; when the shoot created phytomers, their lengths given at once and the shoot's elements updated."""
    initiators = (ShootArchitecture,)

    def __init__(self, data_structure, time_step):
        self.architecture = ShootArchitecture(data_structure=data_structure, time_step=time_step)
        self.development = ScenarioDevelopment(data_structure=data_structure, time_step=time_step)
        self.components = [self.development, self.architecture]

    def run(self):
        ds = self.architecture.data_structure
        self.development()
        before = (ds.topology_version, ds.n_nodes())
        self.architecture()
        if (ds.topology_version, ds.n_nodes()) != before:
            self.development.write()
            self.architecture()


@dataclass
class Planting:
    """A plant of a stand: its position (cm), its azimuth (°) and its delay (°Cd: its scenario starts that much
    later)."""
    x: float = 0.
    y: float = 0.
    azimuth: float = 0.
    delay: float = 0.


def sown_stand(n_plants, n_rows=4, row_spacing=15., plant_spacing=4., max_delay=60., seed=0):
    """A sown stand of *n_plants* (Planting: positions in cm) in *n_rows* rows, a little off their sowing spots, their
    azimuths random, their delays random up to *max_delay* (°Cd)."""
    random = np.random.default_rng(seed)
    per_row = int(np.ceil(n_plants / n_rows))
    return [Planting(x=(k % per_row) * plant_spacing + random.normal(0., 0.5),
                     y=(k // per_row) * row_spacing + random.normal(0., 0.5), azimuth=random.uniform(0., 360.),
                     delay=random.uniform(0., max_delay)) for k in range(n_plants)]


def shoot_scene(plants=None, model=ScenarioShoot, time_step=3600, plant_parameters=None, **parameters):
    """A scene of *plants* ([Planting], by default one at the origin), each a *model* (ScenarioShoot or one extending
    it) with *parameters* (its components' parameters, e.g. axes_table, phytomer_table, leaves, n_sect,
    thermal_time_step, tiller_inclination, leaf_number; a plant's delay is its Planting's) and its own
    *plant_parameters* (one dict per plant, e.g. its leaf_number)."""
    plants = plants or [Planting()]
    own = plant_parameters or [{}] * len(plants)
    table = pd.DataFrame([dict(plant=f"p{k}", model=model, x=plant.x * 0.01, y=plant.y * 0.01, z=0.,
                               rotation=float(plant.azimuth),
                               scenario={"parameters": dict(parameters, delay=float(plant.delay), **own[k])})
                          for k, plant in enumerate(plants)])
    xs, ys = table["x"].to_numpy(), table["y"].to_numpy()
    table.attrs.update(xrange=max(float(np.ptp(xs)), 0.1), yrange=max(float(np.ptp(ys)), 0.1))
    return Scene(table, time_step=time_step)
