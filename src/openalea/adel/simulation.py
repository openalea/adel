"""
Adel's wheat shoots for users: a plant grown in a few lines, its results as tables.

    from openalea.adel.simulation import simulate, parameters

    result = simulate(days=80)                    one plant, Adel's sample plant, at 15 °C
    result.summary                                per day: thermal time, Haun stage, leaves, tillers, height, leaf area
    result.organs()                               per organ at the end: its axis, kind, rank, lengths and area
    result.axes()                                 per axis (main stem, tillers): its phytomers and Haun stage
    result.draw()                                 the plant in 3D
    parameters()                                  the parameters: their value, unit and description
    simulate({"phyllochron": 100.}, days=80)      leaves appearing every 100 °Cd
    simulate(days=80, temperature=series)         a measured air temperature (°C), by day
    simulate(days=80, outputs="outputs")          the run's files too: run.json, scenario.json, summary.csv, the
                                                  plants' summaries (Plant/), a frame of the plants every 10 days
    tables_from_measurements(dimensions)          Adel's tables from a main stem's measured dimensions

The plant grows one step per day; its development follows the thermal time, the sum of the days' temperatures above
*base_temperature* (ScenarioDevelopment's, from the scene's forcing air_temperature). The organs' lengths come from
ScenarioDevelopment (openalea.adel.development), a simple linear scenario, and ShootArchitecture
(openalea.adel.architecture) builds the plant from them. Lengths in m, areas in m2.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from openalea.metafspm.coupling.forcings import as_forcing
from openalea.metafspm.scene.recorder import SceneRecorder
from openalea.metafspm.scene.scenarios import parameter_table

from openalea.adel.architecture import ShootArchitecture
from openalea.adel.development import DAY, ScenarioDevelopment, ScenarioShoot, random_delays, shoot_scene, sown_stand
from openalea.adel.tables import wheat_tables, with_tillers

# The parameters a user changes first: the plant described by four numbers (without tables, instead of Adel's sample
# plant), its development and its senescence
MAIN_PARAMETERS = ("leaf_number", "tiller_number", "longest_blade_length", "longest_blade_width", "phyllochron",
                   "blade_duration", "sheath_duration", "internode_start", "internode_duration", "senescence_start",
                   "senescence_duration", "tiller_inclination")


def parameters(values: dict = None, main_only: bool = False) -> pd.DataFrame:
    """The parameters of the plant model (ShootArchitecture's and ScenarioDevelopment's), one row each: its name, its
    value (the default, or the one in *values*), its unit, its description, its component, and whether it is one of
    the main ones (those first), and its bounds (min_value, max_value). A name both components declare (phyllochron)
    takes one value for both."""
    table = parameter_table(ShootArchitecture, ScenarioDevelopment, values=values, main=MAIN_PARAMETERS,
                            exclude=("thermal_time_step",))      # unused: the temperature gives the degree-days
    return table[table["main"]].reset_index(drop=True) if main_only else table


MEASURED_COLUMNS = {"L_blade": "Ll", "W_blade": "Lw", "L_sheath": "Gl", "W_sheath": "Gd", "L_internode": "El",
                    "W_internode": "Ed"}


def tables_from_measurements(dimensions: pd.DataFrame, tiller_number: int = 4, phyllochron: float = 88.92) -> dict:
    """Adel's tables ({"axes_table", "phytomer_table"}, parameters for simulate) of a plant whose main stem was
    measured: *dimensions* has one row per phytomer from the bottom, its columns plantgen's (index_phytomer, L_blade,
    W_blade, L_sheath, W_sheath, L_internode, W_internode, in cm). The tillers and the leaf shapes are a generated
    plant's (tables.wheat_tables, *tiller_number* tillers): each tiller takes the main stem's ranks, 0.9 times."""
    axes, phytomers = wheat_tables(leaf_number=len(dimensions), tiller_number=tiller_number, phyllochron=phyllochron)
    main = phytomers.loc[phytomers["axe"] == "MS"].copy()
    ordered = dimensions.sort_values("index_phytomer")
    for column, name in MEASURED_COLUMNS.items():
        main[name] = ordered[column].to_numpy(dtype=float)
    return {"axes_table": axes.copy(), "phytomer_table": with_tillers(axes, main)}


class _Summary:
    """The summary row of the day just grown (the scene's iteration - 1), every *every* days and on the last day."""

    def __init__(self, every: int, last_day: int):
        self.every, self.last_day = every, last_day

    def __call__(self, scene):
        day = scene.iteration - 1
        if day < 0 or not (day % self.every == 0 or day == self.last_day):
            return None
        population = scene.populations[0]
        return _summary_row(float(day), population.data_structure, population.instance.architecture)


@dataclass
class Simulation:
    """A grown plant (or stand): its *summary* (a DataFrame, one row per recorded day; with several plants, their
    means) and its *scene* (metafspm's, for more: its DataStructure holds every variable of every element)."""
    scene: object
    summary: pd.DataFrame
    parameters: dict

    @property
    def data_structure(self):
        return self.scene.populations[0].data_structure

    @property
    def architecture(self):
        return self.scene.populations[0].instance.architecture

    def _plant_of(self, scale):
        return self.data_structure.owner("Plant", of=scale)

    def axes(self) -> pd.DataFrame:
        """The axes now, one row each: its plant (from 0), its name (MS the main stem, T1, T2... the tillers), the
        number of phytomers it has, its Haun stage (its leaves appeared) and its final number of leaves."""
        ds = self.data_structure
        return pd.DataFrame(dict(plant=self._plant_of("Axis"), axis=self.architecture.axis_names(),
                                 phytomers=np.asarray(ds.get("axis_phytomers")),
                                 haun_stage=np.asarray(ds.get("haun_stage")),
                                 final_leaf_number=np.asarray(ds.get("final_leaf_number"))))

    def organs(self) -> pd.DataFrame:
        """The organs now, one row each: its plant, its axis, its kind (internode, sheath or blade), its rank on the
        axis (from 1, the bottom), its length, its visible length (out of the sheaths), its senesced length (m) and
        its visible area (m2)."""
        ds = self.data_structure
        return pd.DataFrame(dict(plant=self._plant_of("Organ"),
                                 axis=self.architecture.axis_names()[ds.owner("Axis", of="Organ")],
                                 organ=self.architecture.organ_names(), rank=np.asarray(ds.get("organ_rank")),
                                 length=np.asarray(ds.get("organ_length")),
                                 visible_length=np.asarray(ds.get("visible_length")),
                                 senesced_length=np.asarray(ds.get("organ_senesced_length")),
                                 area=np.asarray(ds.get("organ_area"))))

    def draw(self, ax=None, color_by: str = None, view=(15., -60.)):
        """The plant in 3D: its elements green, or straw when more senesced than green; or coloured by *color_by*
        (a variable of the elements, e.g. "area"). Returns the matplotlib axes."""
        import matplotlib.pyplot as plt
        from openalea.adel.display import state_colors
        from openalea.metafspm.display import draw, value_colors

        if ax is None:
            ax = plt.figure(figsize=(6, 6)).add_subplot(projection="3d")
        ds = self.data_structure
        colors = (state_colors(ds.get("is_green")) if color_by is None
                  else value_colors(np.asarray(ds.get(color_by), dtype=float)))
        draw(ax, self.architecture.triangles(), colors=colors, view=view)
        return ax


def _summary_row(day: float, ds, architecture) -> dict:
    """The plants' means: thermal time, the main stem's Haun stage and visible leaves, tillers, height, areas."""
    n_plants = len(ds.entity_ids("Plant"))
    main_stem = np.asarray(ds.get("bearer_row")) < 0
    blade = architecture.organ_names() == "blade"
    on_main = main_stem[ds.owner("Axis", of="Organ")]
    visible_leaf = blade & on_main & (np.asarray(ds.get("visible_length")) > 0)
    visible_leaves = ds.map_to(visible_leaf.astype(float), "Organ", "Plant", "sum")
    stem = np.asarray(ds.labels()) == "StemElement"
    height = ds.map_to(np.where(stem, np.asarray(ds.get("z2")), 0.), "node", "Plant", "max")    # its stems' top
    leaf = np.asarray(ds.labels()) == "LeafElement"
    return dict(day=day, thermal_time=float(np.mean(ds.get("thermal_time"))),
                haun_stage=float(np.mean(np.asarray(ds.get("haun_stage"))[main_stem])),
                visible_leaves=float(np.mean(visible_leaves)),
                tillers=float(np.mean(np.maximum(np.asarray(ds.get("axis_number")) - 1., 0.))),
                height=float(np.mean(height)),
                green_leaf_area=float(np.sum(np.asarray(ds.get("green_area"))[leaf])) / n_plants,
                senesced_leaf_area=float(np.sum(np.asarray(ds.get("senesced_area"))[leaf])) / n_plants)


def simulate(parameters: dict = None, days: float = 80, temperature=15., base_temperature: float = 0.,
             plants: int = 1, every: int = 1, seed: int = 0, model=ScenarioShoot, outputs: str = None,
             compute=None) -> Simulation:
    """Grow *plants* wheat plants (one by default; several: a stand sown in rows, with random development delays of up
    to 60 °Cd) for *days* days, one step per day, from the plant model's default parameters (Adel's sample plant)
    changed by *parameters* ({name: value}; see parameters()). The summary is recorded every *every* days.

    *temperature* (°C; 15 by default) is the daily mean temperature: a number, a function of the day, or a pandas
    Series indexed by day (the scene's forcing air_temperature). Each day adds its temperature above
    *base_temperature* (°C) to the thermal time (°Cd), which drives the development; the thermal time starts at 0 at
    day 0 (the sowing). *model*: the plant model,
    ScenarioShoot or one extending it (with more components, e.g. a light model).

    *outputs*: a folder for the run's files (metafspm's SceneRecorder): run.json, scenario.json, summary.csv (the
    summary), Plant/<model>/summaries.csv (every variable per plant, every day) and a frame of the plants every 10
    days (frames/, for ParaView or another platform).
    *compute*: the Scene's run configuration (metafspm's Compute: its threads, its memory trim; None: Compute())."""
    values = dict(parameters or {})
    values.setdefault("start", 0.)
    values.setdefault("base_temperature", base_temperature)
    recorder = SceneRecorder(outputs, summaries=outputs is not None,
                             tables={"summary": _Summary(int(every), int(days))},
                             frame_period=10 if outputs is not None else None)
    forcings = {"air_temperature": as_forcing(temperature, time_unit=DAY)}
    if plants == 1:
        built = shoot_scene(model=model, time_step=DAY, forcings=forcings, recorder=recorder, compute=compute,
                            **values)
    else:
        built = shoot_scene(sown_stand(plants, seed=seed), model=model, time_step=DAY,
                            delays=random_delays(plants, seed=seed), forcings=forcings, recorder=recorder,
                            compute=compute, **values)
    built.simulate(int(days) + 1)                      # day 0's step sets the thermal time to its start
    return Simulation(scene=built, summary=recorder.table("summary"), parameters=dict(parameters or {}))
