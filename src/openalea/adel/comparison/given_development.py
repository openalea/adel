"""
The harness' scenario (legacy_run.given_lengths) as a metafspm component: GivenDevelopment writes what a development
and a senescence model write to ShootArchitecture (the coupling contract, openalea.adel.architecture), from a thermal
time advanced at each call, vectorised on the DataStructure:

- per plant, its axis_number: its axes whose first leaf appeared (axes table's emf1);
- per axis, its phytomer_number: its metamers appeared, one per phyllochron from its first leaf, nf at most;
- per organ, its organ_length, visible_length and organ_senesced_length, and a blade's age: a blade grows linearly to
  its final length (phytomer table's Ll) in 1.6 phyllochrons, its sheath (Gl) from 1 phyllochron for 1, its internode
  (El) from 2 phyllochrons for 1.5; the blade and the sheath visible above the previous sheath, the internode half
  visible; the blade senesces from 5 phyllochrons for 3; its age is the thermal time since it appeared.

A plant's thermal time starts at `start` (by default its main stem's first leaf) minus its `delay`. GivenShoot is the
plant model of ShootArchitecture grown by it, for the comparison harness, the examples and the tutorial:

    table = pd.DataFrame([dict(plant="p0", model=GivenShoot, x=0., y=0., z=0., rotation=0.,
                               scenario={"parameters": dict(axes_table=axes, phytomer_table=phytomers, delay=0.)})])
    scene = Scene(table, time_step=3600)
    scene.run()                                        # 10 °Cd (thermal_time_step) of development, then the shoot
"""
from dataclasses import dataclass

import numpy as np

from openalea.adel.architecture import SCALES, ShootArchitecture, doc
from openalea.metafspm.coupling.component import FunctionalComponent, input_variable, parameter, state_variable
from openalea.metafspm.solve.decorator import state

BY = "model_given_development"
ORGAN_KINDS = ("StemInternode", "Sheath", "Blade")


@dataclass
class GivenDevelopment(FunctionalComponent):
    """The harness' organ lengths, given at each call from the plants' thermal time."""

    # --- Read from ShootArchitecture ---
    axis_row: int = input_variable(**doc('adim', "The axis' row in the axes table"), by='model_shoot_architecture',
                                   initialize=-1, scale=SCALES.Axis, dtype='int')
    organ_rank: int = input_variable(**doc('adim', "The organ's phytomer rank"), by='model_shoot_architecture',
                                     initialize=0, scale=SCALES.Organ, dtype='int')

    # --- Written to ShootArchitecture ---
    thermal_time: float = state_variable(**doc('°Cd', "The plant's thermal time in the scenario", "NaN: not started"),
                                         state_variable_type="intensive", initialize=np.nan, scale=SCALES.Plant)
    axis_number: float = state_variable(**doc('adim', "The number of the plant's axes emerged"),
                                        state_variable_type="descriptor", initialize=1., scale=SCALES.Plant)
    phytomer_number: float = state_variable(**doc('adim', "The number of phytomers the axis has"),
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
    axes_table: object = parameter(**doc('', "Adel's axeT (ShootArchitecture's): emf1 and nf of each axis"), by=BY,
                                   default=None)
    phytomer_table: object = parameter(**doc('', "Adel's phytoT (ShootArchitecture's): Ll, Gl, El of each rank"),
                                       by=BY, default=None)
    table_plant: int = parameter(**doc('adim', "The plant's number in the tables"), by=BY, default=1)
    length_unit: float = parameter(**doc('m', "The tables' length unit in m"), by=BY, default=0.01)
    phyllochron: float = parameter(**doc('°Cd', "The phyllochron of every axis"), by=BY, default=88.92)
    thermal_time_step: float = parameter(**doc('°Cd', "The thermal time of a call"), by=BY, default=10.)
    start: float = parameter(**doc('°Cd', "The scenario's start", "NaN: the main stem's first leaf (emf1)"), by=BY,
                             default=np.nan)
    delay: float = parameter(**doc('°Cd', "The plant's delay: its scenario starts that much later"), by=BY,
                             default=0.)

    @state
    def _thermal_time(self, thermal_time) -> tuple[np.ndarray, str, np.ndarray, str, np.ndarray, str, np.ndarray,
                                                   str, np.ndarray, str, np.ndarray, str, np.ndarray]:
        """The plants' thermal time advanced by a step (set at start - delay at the first call), and what is given
        at it."""
        thermal_time = np.asarray(thermal_time, dtype=float)
        first = self._start() - self._plant_values("delay")
        thermal_time = np.where(np.isnan(thermal_time), first, thermal_time + self._plant_values("thermal_time_step"))
        given = self.given(thermal_time)
        return (thermal_time, "axis_number", given["axis_number"], "phytomer_number", given["phytomer_number"],
                "organ_length", given["organ_length"], "visible_length", given["visible_length"],
                "organ_senesced_length", given["organ_senesced_length"], "age", given["age"])

    def write(self):
        """Write what is given at the plants' current thermal time again: the organs created since are given their
        lengths (as cnwgrass writes a new metamer's at once)."""
        ds = self.data_structure
        for name, values in self.given(np.asarray(ds.get("thermal_time"), dtype=float)).items():
            ds.set(name, values)

    def given(self, thermal_time):
        """{variable: values} at the plants' *thermal_time*: their axis numbers, their axes' phytomer numbers, their
        organs' lengths (m) and blades' ages."""
        ds = self.data_structure
        tables = self._tables()
        unit, phyllochron = self._plant_values("length_unit"), self._plant_values("phyllochron")
        table_plant = self._plant_values("table_plant").astype(int)
        emf1, nf = tables["emf1"], tables["nf"]

        # Plants: their axes emerged (in their table's rows)
        axis_number = np.array([float(np.sum(thermal_time[p] >= emf1[tables["plant_rows"][number]]))
                                for p, number in enumerate(table_plant.tolist())])

        # Axes: their metamers appeared
        axis_row = np.asarray(ds.get("axis_row"), dtype=int)
        axis_time = ds.map_to(thermal_time, "Plant", "Axis", "broadcast")
        axis_phyllochron = ds.map_to(phyllochron, "Plant", "Axis", "broadcast")
        ranks = np.arange(1, tables["Ll"].shape[1])
        appeared = emf1[axis_row][:, None] + (ranks[None, :] - 1) * axis_phyllochron[:, None]
        phytomer_number = np.sum((axis_time[:, None] >= appeared) & (ranks[None, :] <= nf[axis_row][:, None]),
                                 axis=1).astype(float)

        # Organs: their lengths
        n_axes = len(ds.entity_ids("Axis"))
        row = axis_row[ds.map_to(np.arange(n_axes, dtype=float), "Axis", "Organ", "broadcast").astype(int)]
        rank = np.asarray(ds.get("organ_rank"), dtype=int)
        time = ds.map_to(thermal_time, "Plant", "Organ", "broadcast")
        organ_phyllochron = ds.map_to(phyllochron, "Plant", "Organ", "broadcast")
        organ_unit = ds.map_to(unit, "Plant", "Organ", "broadcast")
        kind = self._organ_kinds()
        organ_appeared = emf1[row] + (rank - 1) * organ_phyllochron

        def ramp(begin, duration):
            return np.clip((time - begin) / duration, 0., 1.)

        def dimension(name, ranks_of):
            values = tables[name][row, np.clip(ranks_of, 0, None)]
            return np.where((ranks_of >= 1) & np.isfinite(values), values, 0.)

        blade = dimension("Ll", rank) * ramp(organ_appeared, 1.6 * organ_phyllochron)
        sheath = dimension("Gl", rank) * ramp(organ_appeared + organ_phyllochron, organ_phyllochron)
        internode = dimension("El", rank) * ramp(organ_appeared + 2 * organ_phyllochron, 1.5 * organ_phyllochron)
        previous_sheath = dimension("Gl", rank - 1) * ramp(organ_appeared, organ_phyllochron)
        blade_senesced = np.minimum(dimension("Ll", rank) * ramp(organ_appeared + 5 * organ_phyllochron,
                                                                 3 * organ_phyllochron) * (blade > 0), blade)
        length = np.select([kind == 0, kind == 1, kind == 2], [internode, sheath, blade], 0.)
        visible = np.select([kind == 0, kind == 1, kind == 2],
                            [0.5 * internode, np.maximum(0., sheath - previous_sheath),
                             np.maximum(0., blade - previous_sheath)], 0.)
        senesced = np.where(kind == 2, blade_senesced, 0.)
        age = np.where((kind == 2) & (time >= organ_appeared), time - organ_appeared,
                       np.asarray(ds.get("age"), dtype=float))              # a blade not appeared yet keeps its age
        return dict(axis_number=axis_number, phytomer_number=phytomer_number, organ_length=length * organ_unit,
                    visible_length=visible * organ_unit, organ_senesced_length=senesced * organ_unit, age=age)

    def _plant_values(self, name):
        return np.asarray(self.data_structure.parameter_view(name, "Plant"), dtype=float)

    def _tables(self):
        """The axes table's emf1 and nf per row, each table plant's rows, the phytomer table's Ll, Gl, El per (axis
        row, rank) in its unit; read once."""
        cache = self.__dict__.get("_tables_cache")
        if cache is None:
            rows = self.axes_table.reset_index(drop=True)
            dimensions = ShootArchitecture._dimensions(self.phytomer_table, rows, 1.)
            cache = dict(emf1=rows["emf1"].to_numpy(dtype=float), nf=rows["nf"].to_numpy(dtype=int),
                         plant_rows={int(number): ShootArchitecture._plant_rows(rows, int(number))
                                     for number in set(rows["plant"])},
                         Ll=dimensions["Ll"], Gl=dimensions["Gl"], El=dimensions["El"])
            self.__dict__["_tables_cache"] = cache
        return cache

    def _start(self):
        """Each plant's start: its parameter, by default its main stem's first leaf (emf1)."""
        tables = self._tables()
        start = self._plant_values("start")
        main_stem = np.array([tables["emf1"][tables["plant_rows"][number][0]]
                              for number in self._plant_values("table_plant").astype(int).tolist()])
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


class GivenShoot:
    """A plant model: ShootArchitecture grown by GivenDevelopment. At each call the development first, then the shoot;
    when the shoot created phytomers, their lengths given at once and the shoot's elements updated."""
    initiators = (ShootArchitecture,)

    def __init__(self, data_structure, time_step):
        self.architecture = ShootArchitecture(data_structure=data_structure, time_step=time_step)
        self.development = GivenDevelopment(data_structure=data_structure, time_step=time_step)
        self.components = [self.development, self.architecture]

    def run(self):
        ds = self.architecture.data_structure
        self.development()
        before = (ds.topology_version, ds.n_nodes())
        self.architecture()
        if (ds.topology_version, ds.n_nodes()) != before:
            self.development.write()
            self.architecture()
