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
    tables_from_measurements(dimensions)          Adel's tables from a main stem's measured dimensions

The plant grows one step per day; its development follows the thermal time, the sum of the days' temperatures above
*base_temperature*. The organs' lengths come from ScenarioDevelopment (openalea.adel.development), a simple linear
scenario, and ShootArchitecture (openalea.adel.architecture) builds the plant from them. Lengths in m, areas in m2.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from openalea.adel.architecture import ShootArchitecture
from openalea.adel.development import ScenarioDevelopment, ScenarioShoot, random_delays, shoot_scene, sown_stand
from openalea.adel.tables import wheat_tables, with_tillers

DAY = 86400.

# The parameters a user changes first: the plant described by four numbers (without tables, instead of Adel's sample
# plant), its development and its senescence
MAIN_PARAMETERS = ("leaf_number", "tiller_number", "longest_blade_length", "longest_blade_width", "phyllochron",
                   "blade_duration", "sheath_duration", "internode_start", "internode_duration", "senescence_start",
                   "senescence_duration", "tiller_inclination")


def parameters(values: dict = None, main_only: bool = False) -> pd.DataFrame:
    """The parameters of the plant model (ShootArchitecture's and ScenarioDevelopment's), one row each: its name, its
    value (the default, or the one in *values*), its unit, its description, its component, and whether it is one of
    the main ones (those first). A name both components declare (phyllochron) takes one value for both."""
    tables = []
    for component in (ShootArchitecture, ScenarioDevelopment):
        declarations = component.declarations()
        table = declarations[declarations["kind"].isin(("parameter", "setting"))][["name", "default", "unit",
                                                                                     "description"]]
        tables.append(table.assign(component=component.__name__))
    table = pd.concat(tables).drop_duplicates("name", keep="last").rename(columns={"default": "value"})
    table = table[table["name"] != "thermal_time_step"].reset_index(drop=True)      # simulate gives it each day
    for name, value in (values or {}).items():
        table.loc[table["name"] == name, "value"] = value
    table["main"] = table["name"].isin(MAIN_PARAMETERS)
    order = {name: rank for rank, name in enumerate(MAIN_PARAMETERS)}
    table = table.assign(rank=table["name"].map(order).fillna(len(order))).sort_values("rank", kind="stable")
    table = table.drop(columns="rank").reset_index(drop=True)
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


def _value_at(value, day: float) -> float:
    """A given input at *day*: a number; a function of the day; a pandas Series indexed by day, or a sequence (one
    value per day from day 0), linearly interpolated."""
    if callable(value):
        return float(value(day))
    if isinstance(value, pd.Series):
        return float(np.interp(day, value.index.to_numpy(dtype=float), value.to_numpy(dtype=float)))
    if np.ndim(value) > 0:
        values = np.asarray(value, dtype=float)
        return float(np.interp(day, np.arange(values.size), values))
    return float(value)


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
        ds = self.data_structure
        return ds.map_to(np.arange(len(ds.entity_ids("Plant")), dtype=float), "Plant", scale, "broadcast").astype(int)

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
        axis_of = ds.map_to(np.arange(len(ds.entity_ids("Axis")), dtype=float), "Axis", "Organ", "broadcast")
        return pd.DataFrame(dict(plant=self._plant_of("Organ"),
                                 axis=self.architecture.axis_names()[axis_of.astype(int)],
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
    on_main = ds.map_to(main_stem.astype(float), "Axis", "Organ", "broadcast") > 0
    plant_of_organ = ds.map_to(np.arange(n_plants, dtype=float), "Plant", "Organ", "broadcast").astype(int)
    visible_leaves = np.bincount(plant_of_organ, weights=(blade & on_main
                                                           & (np.asarray(ds.get("visible_length")) > 0)),
                                 minlength=n_plants)
    stem = np.isin(np.asarray(ds.labels()), ("StemElement",))
    plant_of_node = ds.map_to(np.arange(n_plants, dtype=float), "Plant", "node", "broadcast").astype(int)
    top = np.where(stem, np.asarray(ds.get("z2")), -np.inf)
    height = np.array([max(0., top[plant_of_node == p].max(initial=0.)) for p in range(n_plants)])
    leaf = np.asarray(ds.labels()) == "LeafElement"
    return dict(day=day, thermal_time=float(np.mean(ds.get("thermal_time"))),
                haun_stage=float(np.mean(np.asarray(ds.get("haun_stage"))[main_stem])),
                visible_leaves=float(np.mean(visible_leaves)),
                tillers=float(np.mean(np.maximum(np.asarray(ds.get("axis_number")) - 1., 0.))),
                height=float(np.mean(height)),
                green_leaf_area=float(np.sum(np.asarray(ds.get("green_area"))[leaf])) / n_plants,
                senesced_leaf_area=float(np.sum(np.asarray(ds.get("senesced_area"))[leaf])) / n_plants)


def simulate(parameters: dict = None, days: float = 80, temperature=15., base_temperature: float = 0.,
             plants: int = 1, every: int = 1, seed: int = 0, model=ScenarioShoot) -> Simulation:
    """Grow *plants* wheat plants (one by default; several: a stand sown in rows, with random development delays of up
    to 60 °Cd) for *days* days, one step per day, from the plant model's default parameters (Adel's sample plant)
    changed by *parameters* ({name: value}; see parameters()). The summary is recorded every *every* days.

    *temperature* (°C; 15 by default) is the daily mean temperature: a number, a function of the day, or a pandas
    Series indexed by day. Each day adds its temperature above *base_temperature* (°C) to the thermal time (°Cd),
    which drives the development; the thermal time starts at 0 at day 0 (the sowing). *model*: the plant model,
    ScenarioShoot or one extending it (with more components, e.g. a light model)."""
    values = dict(parameters or {})
    values.setdefault("start", 0.)
    if plants == 1:
        built = shoot_scene(model=model, time_step=DAY, **values)
    else:
        built = shoot_scene(sown_stand(plants, seed=seed), model=model, time_step=DAY,
                            delays=random_delays(plants, seed=seed), **values)
    ds, architecture = built.populations[0].data_structure, built.populations[0].instance.architecture
    n_plants = len(ds.entity_ids("Plant"))
    built.run()                                        # day 0: the thermal time set to its start
    rows = [_summary_row(0., ds, architecture)]
    for day in range(1, int(days) + 1):
        degrees = max(0., _value_at(temperature, day - 1) - base_temperature)
        ds.set("thermal_time_step", np.full(n_plants, degrees))
        built.run()
        if day % every == 0 or day == int(days):
            rows.append(_summary_row(float(day), ds, architecture))
    return Simulation(scene=built, summary=pd.DataFrame(rows), parameters=dict(parameters or {}))
