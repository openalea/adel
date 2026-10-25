"""
The tables ShootArchitecture reads: Adel's axes table (axeT: one row per axis of each plant) and phytomer table
(phytoT: one row per axis and rank), lengths in the tables' unit (Adel's cm). A plant's tables are, in order of
precedence:

- the tables given (ShootArchitecture's axes_table and phytomer_table, the plant's rows by its table_plant);
- tables generated from a few scalars, when any is given (wheat_tables: ShootArchitecture's per-plant parameters
  leaf_number, tiller_number, longest_blade_length, longest_blade_width, with its phyllochron);
- Adel's sample plant (sample_tables).

    axes, phytomers = sample_tables()
    axes, phytomers = wheat_tables(leaf_number=11, tiller_number=2, longest_blade_length=0.2,
                                   longest_blade_width=0.015, phyllochron=100.)

Each plant's tables are kept on the MPG at its initiation (record_plant_tables); population_tables gathers them into
one pair of tables for the population, every plant's rows at an offset: the rows the components index (axis_row).
"""
import os
from functools import lru_cache

import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(__file__), "data")
PROFILES = os.path.join(DATA, "Wheat_dimension_profiles.csv")
TILLER_FACTOR = 0.9                                    # a tiller's dimensions: its bearer's main stem's, times it
SCALARS = ("leaf_number", "tiller_number", "longest_blade_length", "longest_blade_width")
SCALAR_DEFAULTS = dict(leaf_number=12, tiller_number=4, longest_blade_length=0.19, longest_blade_width=0.016)


def sample_tables(tiller_factor=TILLER_FACTOR):
    """Adel's sample plant (data/axeTable.csv: a main stem and four tillers; data/phytoT.csv: the main stem's ranks),
    each tiller given the main stem's first nf ranks, its lengths and diameters *tiller_factor* times."""
    axes = pd.read_csv(os.path.join(DATA, "axeTable.csv"))
    phytomers = pd.read_csv(os.path.join(DATA, "phytoT.csv"))
    return axes, with_tillers(axes, phytomers, tiller_factor)


def with_tillers(axes, phytomers, tiller_factor=TILLER_FACTOR):
    """*phytomers* with the ranks of each tiller of *axes* it lacks: its plant's main stem's first nf ranks, lengths
    and diameters *tiller_factor* times."""
    rows = [phytomers]
    for _, axis in axes.loc[axes["axe"] != "MS"].iterrows():
        present = (phytomers["plant"] == axis["plant"]) & (phytomers["axe"] == axis["axe"])
        if present.any():
            continue
        main = phytomers.loc[(phytomers["plant"] == axis["plant"]) & (phytomers["axe"] == "MS")]
        tiller = main.loc[main["n"] <= int(axis["nf"])].copy()
        tiller["axe"] = axis["axe"]
        for name in ("Ll", "Lw", "Gl", "Gd", "El", "Ed"):
            tiller[name] = tiller[name] * tiller_factor
        rows.append(tiller)
    return pd.concat(rows, ignore_index=True)


@lru_cache(maxsize=None)
def _profiles():
    return pd.read_csv(PROFILES)


@lru_cache(maxsize=64)
def wheat_tables(leaf_number=12, tiller_number=4, longest_blade_length=0.19, longest_blade_width=0.016,
                 phyllochron=88.92, unit=0.01, shape_classes=6, dynamic_shapes=False):
    """(axes table, phytomer table) of a wheat plant from a few scalars (lengths in m, the tables in *unit*):

    - the main stem: *leaf_number* leaves, its first at 0 °Cd; its blades', sheaths' and internodes' lengths and
      widths along Adel's reference profiles (data/Wheat_dimension_profiles.csv, plantgen's WheatDimensions: by the
      rank relative to the longest leaf's, 2.4 below the flag leaf's), scaled so that its longest blade is
      *longest_blade_length* long and *longest_blade_width* wide;
    - *tiller_number* tillers (T1 on the main stem's phytomer 1, T2 on its 2nd, ...): T_k appears k + 0.6
      phyllochrons after the main stem, with leaf_number - 2 - k leaves, its ranks the main stem's times 0.9; inclined
      60° at their base, straightening 4 cm from their insertion (used with ShootArchitecture's tiller_inclination);
    - the shape classes (Lindex) by the leaf's rank from the top: the samples' (class min(ntop, *shape_classes*)), or
      for a dynamic database (*dynamic_shapes*, echap's) echap's (echap_shape_indices); shape 1 in its class
      (Lseed).

    Cached: the same scalars give the same tables (the same objects)."""
    leaf_number, tiller_number = int(leaf_number), int(tiller_number)
    profiles = _profiles()
    ranks = np.arange(1, leaf_number + 1)
    relative = np.clip(ranks - (leaf_number - 2.4), profiles["xn"].min(), profiles["xn"].max())

    def profile(column):
        return np.interp(relative, profiles["xn"], profiles[column])

    blade_length = profile("L_blade")
    length_scale = longest_blade_length / unit / blade_length.max()
    width_scale = longest_blade_width / unit / profile("W_blade").max()
    ntop = leaf_number - ranks + 1
    shape_class = echap_classes(ranks, leaf_number) if dynamic_shapes else np.minimum(ntop, shape_classes)
    main = pd.DataFrame(dict(plant=1, axe="MS", n=ranks, Ll=blade_length * length_scale,
                             Lw=profile("W_blade") * width_scale, Gl=profile("L_sheath") * length_scale,
                             Gd=profile("W_sheath") * width_scale, El=profile("L_internode") * length_scale,
                             Ed=profile("W_internode") * width_scale, Azim=180., Lindex=shape_class.astype(float),
                             Lseed=1.))
    axes = [dict(plant=1, axe="MS", nf=float(leaf_number), emf1=0., HS_final=float(leaf_number), hasEar=1, incT=0.,
                 dredT=0., azT=0., azTb=0.)]
    for k in range(1, tiller_number + 1):
        tiller_leaves = leaf_number - 2 - k
        if tiller_leaves < 1:
            break
        axes.append(dict(plant=1, axe=f"T{k}", nf=float(tiller_leaves), emf1=(k + 0.6) * phyllochron,
                         HS_final=float(tiller_leaves), hasEar=0, incT=60., dredT=0.04 / unit, azT=75., azTb=0.))
    axes = pd.DataFrame(axes)
    return axes, with_tillers(axes, main)


def echap_classes(ranks, leaf_number, top_leaves=4):
    """echap's shape classes (echap_leaf.geoLeaf): 2 for the top leaves (leaf_number - rank <= *top_leaves*), else
    1."""
    return np.where(leaf_number - np.asarray(ranks) <= top_leaves, 2, 1)


def echap_shape_indices(axes, phytomers):
    """*phytomers* with echap's shape indices (echap_classes by the plant's main stem's leaf number; shape 1)."""
    adapted = phytomers.copy()
    main_stem = axes.loc[axes["axe"] == "MS"].set_index("plant")["nf"]
    adapted["Lindex"] = echap_classes(adapted["n"], adapted["plant"].map(main_stem).to_numpy(dtype=float))
    adapted["Lseed"] = 1
    return adapted


def plant_tables(parameters, defaults):
    """(axes table, phytomer table, table plant) of a plant from its scenario *parameters* (ShootArchitecture's; its
    class' values in *defaults*): the tables given, else generated when any of its scalars is given, else the
    sample plant's (with echap's shape indices for a dynamic database)."""
    if parameters.get("axes_table") is not None:
        return (parameters["axes_table"], parameters["phytomer_table"],
                int(parameters.get("table_plant", defaults["table_plant"])))
    given = {name: parameters.get(name, np.nan) for name in SCALARS}
    leaves = parameters.get("leaves")
    if all(pd.isna(value) for value in given.values()):
        return (*_sample_tables_once(bool(getattr(leaves, "dynamic", False))), 1)
    values = {name: SCALAR_DEFAULTS[name] if pd.isna(value) else float(value) for name, value in given.items()}
    axes, phytomers = wheat_tables(int(values["leaf_number"]), int(values["tiller_number"]),
                                   values["longest_blade_length"], values["longest_blade_width"],
                                   float(parameters.get("phyllochron", defaults["phyllochron"])),
                                   float(parameters.get("length_unit", defaults["length_unit"])),
                                   len(leaves.classes) if leaves is not None else 6,
                                   bool(getattr(leaves, "dynamic", False)))
    return axes, phytomers, 1


@lru_cache(maxsize=2)
def _sample_tables_once(dynamic_shapes):
    axes, phytomers = sample_tables()
    return axes, echap_shape_indices(axes, phytomers) if dynamic_shapes else phytomers


def record_plant_tables(g, plant, tables):
    """Keep a plant's (axes table, phytomer table, table plant) on the MPG *g*."""
    g.__dict__.setdefault("_shoot_tables", {})[plant] = tables


def population_tables(g, plants):
    """The population's tables from its plants' (record_plant_tables): (axes table, phytomer table, per plant its
    number in them, per plant its rows' offset). The tables of several plants are gathered once (the same objects);
    each source's plants are renumbered source index * 1000 + their number. Cached on *g* for the same plants."""
    plants = tuple(int(plant) for plant in plants)
    cache = g.__dict__.get("_shoot_population_tables")
    if cache is not None and cache[0] == plants:
        return cache[1]
    recorded = g.__dict__.get("_shoot_tables", {})
    sources, numbers, offsets, axes_parts, phytomer_parts = {}, [], [], [], []
    row = 0
    for plant in plants:
        axes, phytomers, number = recorded[plant]
        key = (id(axes), id(phytomers))
        if key not in sources:
            index = len(sources)
            sources[key] = (index, row)
            axes_parts.append(axes.reset_index(drop=True).assign(plant=axes["plant"].to_numpy() + index * 1000))
            phytomer_parts.append(phytomers.assign(plant=phytomers["plant"].to_numpy() + index * 1000))
            row += len(axes)
        index, offset = sources[key]
        numbers.append(index * 1000 + number)
        offsets.append(offset)
    tables = (pd.concat(axes_parts, ignore_index=True), pd.concat(phytomer_parts, ignore_index=True),
              np.array(numbers), np.array(offsets))
    g.__dict__["_shoot_population_tables"] = (plants, tables)
    return tables
