"""
ShootArchitecture (openalea.adel.architecture) driven on the legacy harness' scenario (legacy_run.given_lengths): the
phytomer numbers and organ lengths given at each step, as cnwgrass' development and senescence would write them.

    run = port_run(n_steps=40, leaves="echap")
    run.scene            the metafspm Scene (one plant: ShootArchitecture alone)
    run.steps[i]         a structure record after step run.recorded_steps[i] (structure())

The tables are the legacy harness' (harness_tables).
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from openalea.adel.architecture import ShootArchitecture
from openalea.adel.leaf_shapes import LeafShapes
from openalea.adel.comparison.legacy_run import PHYLLOCHRON, given_lengths, harness_tables, start
from openalea.metafspm.scene.scene import Scene

ORGAN_NAMES = {"StemInternode": "internode", "Sheath": "sheath", "Blade": "blade"}
ELEMENT_VARIABLES = ("length", "area", "green_length", "green_area", "senesced_length", "senesced_area", "is_green",
                     "rolled_part", "rolled_diameter", "width", "flat_base", "flat_top")
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


class ShootOnly:
    """A plant model of ShootArchitecture alone, its inputs given."""
    initiators = (ShootArchitecture,)

    def __init__(self, data_structure, time_step):
        self.architecture = ShootArchitecture(data_structure=data_structure, time_step=time_step)
        self.components = [self.architecture]

    def run(self):
        self.architecture()


@dataclass
class PortRun:
    scene: Scene
    steps: list = field(default_factory=list)
    recorded_steps: list = field(default_factory=list)


def leaf_shapes(leaves):
    return LeafShapes.echap() if leaves == "echap" else LeafShapes.samples()


def port_scene(leaves="samples", time_step=3600):
    axes, phytomers = harness_tables(leaves)
    table = pd.DataFrame([dict(plant="p0", model=ShootOnly, x=0., y=0., z=0., rotation=0.,
                               scenario={"parameters": dict(axes_table=axes, phytomer_table=phytomers,
                                                            table_plant=1, leaves=leaf_shapes(leaves))})])
    table.attrs.update(xrange=1., yrange=1.)
    return Scene(table, time_step=time_step)


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


def write_given(scene, lengths, unit=0.01):
    """Each axis' phytomer number, and its organs' given lengths (cm, legacy_run.given_lengths) in m."""
    ds = scene.populations[0].data_structure
    axes = scene.populations[0].instance.architecture.axes_table
    ds.set("axis_number", np.full(len(ds.entity_ids("Plant")), float(sum(len(lengths[name]) > 0
                                                                         for name in axes["axe"]))))
    names = axis_names(scene)
    ds.set("phytomer_number", np.array([len(lengths[name]) for name in names], dtype=float))
    organ_axis = ds.map_to(np.arange(names.size, dtype=float), "Axis", "Organ", "broadcast").astype(int)
    ranks, organs = np.asarray(ds.get("organ_rank")), organ_names(ds)
    for variable, given in (("organ_length", "length"), ("visible_length", "visible_length"),
                            ("organ_senesced_length", "senesced_length"), ("age", "age")):
        values = np.array(ds.get(variable), dtype=float)
        for i in range(values.size):
            organ = lengths[names[organ_axis[i]]].get(int(ranks[i]), {}).get(organs[i])
            if organ is not None and given in organ:
                values[i] = organ[given] * (unit if given != "age" else 1.)
        ds.set(variable, values)


def structure(scene):
    """{(axis, rank, organ): {"labels": [element labels], "elements": {label: {variable: value}}, variable:
    value}}; leaf elements labelled as the legacy's, by their index (LeafElement1, ...)."""
    ds = scene.populations[0].data_structure
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
    for i, organ in enumerate(organ_ids.tolist()):
        elements = [label_name[g.label(e)] for e in g.components(organ)]
        values = {}
        for e in g.components(organ):
            j = node_index[e]
            label = label_name[g.label(e)]
            if label == "LeafElement":
                label += str(int(element_arrays["element_index"][j]))
            values[label] = {name: element_arrays[name][j].item() for name in ELEMENT_VARIABLES}
        out[(organ_axis[i], int(ranks[i]), names[i])] = dict(labels=elements, elements=values,
                                                              **{name: arrays[name][i].item() for name in variables},
                                                              organ_area=float(np.asarray(ds.get("organ_area"))[i]))
    return out


def port_run(n_steps=40, thermal_time_step=10., leaves="samples", record_every=1, phyllochron=PHYLLOCHRON):
    scene = port_scene(leaves)
    axes, phytomers = harness_tables(leaves)
    run = PortRun(scene=scene)
    for step in range(n_steps + 1):
        lengths = given_lengths(start(axes) + step * thermal_time_step, axes, phytomers, phyllochron=phyllochron)
        write_given(scene, lengths)
        scene.run()
        write_given(scene, lengths)                    # the organs created in the step: their lengths at once,
        run_stage(scene, "postsegmentation")          # as the legacy sets a new metamer's (the elements again)
        if step % record_every == 0 or step == n_steps:
            run.steps.append(structure(scene))
            run.recorded_steps.append(step)
    return run
