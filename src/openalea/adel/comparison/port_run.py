"""
ShootArchitecture (openalea.adel.architecture) driven on the legacy harness' scenario (legacy_run.given_lengths): the
phytomer numbers and organ lengths given at each step by GivenDevelopment (given_development), as cnwgrass' development
and senescence would write them.

    run = port_run(n_steps=40, leaves="echap")
    run.scene            the metafspm Scene (GivenShoot plants: GivenDevelopment and ShootArchitecture)
    run.steps[i]         a structure record after step run.recorded_steps[i] (structure())

The tables are the legacy harness' (harness_tables).
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.comparison.given_development import GivenShoot
from openalea.adel.comparison.legacy_run import PHYLLOCHRON, Planting, harness_tables
from openalea.metafspm.scene.scene import Scene

ORGAN_NAMES = {"StemInternode": "internode", "Sheath": "sheath", "Blade": "blade"}
ELEMENT_VARIABLES = ("length", "area", "green_length", "green_area", "senesced_length", "senesced_area", "is_green",
                     "rolled_part", "rolled_diameter", "width", "flat_base", "flat_top", "x1", "y1", "z1", "x2", "y2",
                     "z2")
STAGE_INDEX = {"potential": 1, "actual": 4, "segmentation": 5, "postsegmentation": 6}


def run_stage(scene, stage):
    """Run only ShootArchitecture's steps of one stage (its row in the Choregrapher's growth schedule)."""
    architecture = scene.populations[0].instance.architecture
    architecture.pull_available_inputs()
    _, groups = architecture.__dict__["_choregraphy"]
    for priority, functors in groups.items():
        if eval(priority)[3] == STAGE_INDEX[stage]:
            for functor in functors:
                functor()


@dataclass
class PortRun:
    scene: Scene
    steps: list = field(default_factory=list)
    recorded_steps: list = field(default_factory=list)


def sown_stand(n_plants, n_rows=4, row_spacing=15., plant_spacing=4., max_delay=60., seed=0):
    """A sown stand of *n_plants* (Planting: positions in cm) in *n_rows* rows, a little off their sowing spots, their
    azimuths random, their delays random up to *max_delay* (°Cd)."""
    random = np.random.default_rng(seed)
    per_row = int(np.ceil(n_plants / n_rows))
    return [Planting(x=(k % per_row) * plant_spacing + random.normal(0., 0.5),
                     y=(k // per_row) * row_spacing + random.normal(0., 0.5), azimuth=random.uniform(0., 360.),
                     delay=random.uniform(0., max_delay)) for k in range(n_plants)]


def leaf_shapes(leaves):
    """The legacy harness' databases (legacy_run.leaf_database): the samples, or echap's with Soissons' midribs."""
    return LeafShapes.echap(xy_model="Soissons_byleafclass") if leaves == "echap" else LeafShapes.samples()


def port_scene(leaves="samples", time_step=3600, rotation=0., phyllochron=PHYLLOCHRON, n_sect=1, plants=None,
               thermal_time_step=10., model=GivenShoot):
    """The scene of the harness' plants (*plants*: [Planting], positions in cm, delays in °Cd; by default one plant at
    the origin, its rotation *rotation*), every plant on the tables' plant 1, the leaf database fitted once for all;
    each scene.run() is *thermal_time_step* °Cd of the scenario (the first one its start). *model*: the plant model,
    GivenShoot or one extending it."""
    axes, phytomers = harness_tables(leaves)
    plants = plants or [Planting(azimuth=rotation)]
    shapes = leaf_shapes(leaves)
    parameters = dict(axes_table=axes, phytomer_table=phytomers, table_plant=1, leaves=shapes, phyllochron=phyllochron,
                      tiller_phyllochron=phyllochron, n_sect=n_sect, thermal_time_step=thermal_time_step)
    table = pd.DataFrame([dict(plant=f"p{k}", model=model, x=plant.x * 0.01, y=plant.y * 0.01, z=0.,
                               rotation=float(plant.azimuth),
                               scenario={"parameters": dict(parameters, delay=float(plant.delay))})
                          for k, plant in enumerate(plants)])
    table.attrs.update(xrange=1., yrange=1.)
    return Scene(table, time_step=time_step)


def plant_of(ds, scale):
    """Each entity's plant (its index among the plants) at *scale*."""
    plants = np.arange(len(ds.entity_ids("Plant")), dtype=float)
    return ds.map_to(plants, "Plant", scale, "broadcast").astype(int)


def organ_names(ds):
    """Each Organ's name (internode, sheath, blade) from its MPG label."""
    g = ds._mtg
    by_value = {getattr(g.labels.Organ, name): organ for name, organ in ORGAN_NAMES.items()}
    return np.array([by_value.get(g.label(int(v))) for v in ds.entity_ids("Organ")], dtype=object)


def axis_names(scene):
    """Each Axis' name (MS, T1, ...) from its row in the axes table."""
    ds = scene.populations[0].data_structure
    axes = scene.populations[0].instance.architecture.axes_table.reset_index(drop=True)
    return axes["axe"].to_numpy()[np.asarray(ds.get("axis_row"))]


def write_inclinations(scene, inclinations):
    """The internodes' inclinations, {(axis, rank): degrees}."""
    ds = scene.populations[0].data_structure
    names = axis_names(scene)
    organ_axis = names[ds.map_to(np.arange(names.size, dtype=float), "Axis", "Organ", "broadcast").astype(int)]
    ranks, organs = np.asarray(ds.get("organ_rank")), organ_names(ds)
    values = np.array(ds.get("inclination"), dtype=float)
    for i in range(values.size):
        if organs[i] == "internode":
            values[i] = inclinations.get((organ_axis[i], int(ranks[i])), values[i])
    ds.set("inclination", values)


def structure(scene, plant=0):
    """{(axis, rank, organ): {"labels": [element labels], "elements": {label: {variable: value, "triangles": (T, 3,
    3)}}, variable: value}} of the plant *plant* (its index); leaf elements labelled as the legacy's, by their index
    (LeafElement1, ...)."""
    ds = scene.populations[0].data_structure
    triangles = scene.populations[0].instance.architecture.triangles()
    g = ds._mtg
    label_name = {}
    for name in ("HiddenElement", "StemElement", "LeafElement"):
        label_name[getattr(g.labels.SubOrgan, name)] = name
    out = {}
    names, ranks = organ_names(ds), np.asarray(ds.get("organ_rank"))
    organ_axis = axis_names(scene)[ds.map_to(np.arange(len(ds.entity_ids("Axis")), dtype=float), "Axis", "Organ",
                                             "broadcast").astype(int)]
    organ_ids = np.asarray(ds.entity_ids("Organ"))
    variables = ("diameter", "shape_mature_length", "shape_max_width", "shape_class", "shape_seed", "ntop",
                 "organ_azimuth")
    arrays = {name: np.asarray(ds.get(name)) for name in variables}
    node_index = dict(zip(np.asarray(ds.entity_ids("node")).tolist(), range(ds.n_nodes())))
    element_arrays = {name: np.asarray(ds.get(name)) for name in ELEMENT_VARIABLES + ("element_index",)}
    organ_plant = plant_of(ds, "Organ")
    for i, organ in enumerate(organ_ids.tolist()):
        if organ_plant[i] != plant:
            continue
        elements = [label_name[g.label(e)] for e in g.components(organ)]
        values = {}
        for e in g.components(organ):
            j = node_index[e]
            label = label_name[g.label(e)]
            if label == "LeafElement":
                label += str(int(element_arrays["element_index"][j]))
            values[label] = {name: element_arrays[name][j].item() for name in ELEMENT_VARIABLES}
            values[label]["triangles"] = triangles.vertices[triangles.owner == j]
        out[(organ_axis[i], int(ranks[i]), names[i])] = dict(labels=elements, elements=values,
                                                              **{name: arrays[name][i].item() for name in variables},
                                                              organ_area=float(np.asarray(ds.get("organ_area"))[i]))
    return out


def port_run(n_steps=40, thermal_time_step=10., leaves="samples", record_every=1, phyllochron=PHYLLOCHRON,
             plant_azimuth=0., inclinations=None, n_sect=1, plants=None):
    """The legacy harness' run (legacy_run's options) on ShootArchitecture; with *plants*, run.steps[i] is a list of
    the plants' structure records."""
    plants_given = plants
    plants = plants or [Planting(azimuth=plant_azimuth)]
    scene = port_scene(leaves, phyllochron=phyllochron, n_sect=n_sect, plants=plants,
                       thermal_time_step=thermal_time_step)
    run = PortRun(scene=scene)
    for step in range(n_steps + 1):
        scene.run()
        if inclinations:
            write_inclinations(scene, inclinations)
            run_stage(scene, "postsegmentation")      # as the legacy sets a new metamer's (the elements again)
        if step % record_every == 0 or step == n_steps:
            records = [structure(scene, plant) for plant in range(len(plants))]
            run.steps.append(records if plants_given else records[0])
            run.recorded_steps.append(step)
    return run
