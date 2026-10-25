"""
The legacy's dynamic geometry (openalea.adel.legacy.adel_dynamic.AdelDyn) driven as cnwgrass drives it, on a
deterministic scenario of organ lengths given at each step (the role of cnwgrass' morphogenesis and senescence), and
its records: the reference ShootArchitecture is compared with.

    run = legacy_run(n_steps=40, leaves="echap")
    run.adel             the legacy AdelDyn
    run.g                its MTG after the last step
    run.steps[i]         a Record after step run.recorded_steps[i]: every vertex of the MTG (scales 1 to 5), its
                         label, scale, complex, topological parent and edge type, its properties (geometry left out),
                         and each element's triangles in the scene's frame (PlantGL's tesselation of its geometry)

The tables (harness_tables): Adel's sample plant (tables.sample_tables: the main stem and 4 tillers, each tiller given
the main stem's first nf ranks, its lengths and diameters 0.9 times).

The scenario (given_lengths), from the main stem's first leaf (its emf1), on every axis: metamer n appears at the axis' first leaf (axeTable's emf1) plus (n - 1)
phyllochrons (add_metamer; a tiller's axis added with its first metamer, add_axe); its blade grows linearly to its
final length (phytoT's Ll) in 1.6 phyllochrons, its sheath (Gl) from 1 phyllochron for 1, its internode (El) from 2
phyllochrons for 1.5; the blade visible above the previous sheath, the sheath above it too, the internode half
visible; the blade senesces from 5 phyllochrons for 3; its age is the thermal time since it appeared. Lengths in the
legacy's unit (cm, phytoT's), thermal time in °Cd.
"""
from dataclasses import dataclass, field

import numpy as np

import openalea.plantgl.all as pgl
from openalea.adel.development import Planting
from openalea.adel.legacy import data_samples
from openalea.adel.tables import echap_shape_indices, sample_tables
from openalea.adel.legacy.adel_dynamic import AdelDyn
from openalea.adel.legacy.echap_leaf import echap_leaves
from openalea.adel.legacy.mtg_editions import add_axe, find_metamers
import pandas as pd

PHYLLOCHRON = 88.92                                    # °Cd, the legacy's default for the main stem
T1_OFFSET = 0.0285                                     # cm, AdelDyn.update_geometry's: T1's elements raised by it


@dataclass
class Record:
    vertices: list
    label: dict
    scale: dict
    complex: dict
    parent: dict
    edge_type: dict
    props: dict
    triangles: dict                                    # {element vid: (T, 3, 3) array}


@dataclass
class LegacyRun:
    adel: AdelDyn
    g: object = None
    steps: list = field(default_factory=list)
    recorded_steps: list = field(default_factory=list)


def leaf_database(name):
    """The leaf shapes: "samples" (data_samples.leaves(), static) or "echap" (echap_leaves, the shapes by age class
    cnwgrass uses)."""
    if name == "samples":
        return data_samples.leaves()
    if name == "echap":
        return echap_leaves(xy_model="Soissons_byleafclass")
    raise ValueError(f"unknown leaf database '{name}'")


def harness_tables(leaves="samples", tiller_factor=0.9):
    """(axes table, phytomer table): Adel's sample plant (tables.sample_tables: each tiller with the main stem's first
    nf ranks, its lengths and diameters *tiller_factor* times), the shape indices those of the leaf database."""
    axes, phytomers = sample_tables(tiller_factor)
    return axes, echap_shape_indices(axes, phytomers) if leaves == "echap" else phytomers


def start(axes):
    """The scenario's start: the main stem's first leaf (its emf1, °Cd)."""
    return float(axes.loc[axes["axe"] == "MS", "emf1"].iloc[0])


def given_lengths(thermal_time, axes, phytoT, phyllochron=PHYLLOCHRON):
    """{axis: {rank n: {organ: {property: value}}}} at *thermal_time*: the metamers appeared on each axis."""
    return {axis["axe"]: axis_lengths(thermal_time, int(axis["nf"]), float(axis["emf1"]), phytoT, axis["axe"],
                                      phyllochron) for _, axis in axes.iterrows()}


def axis_lengths(thermal_time, n_metamers, first_leaf, phytoT, axe, phyllochron=PHYLLOCHRON):
    """{rank n: {organ: {property: value}}} of one axis at *thermal_time* (rank n from 1)."""
    dims = phytoT.loc[phytoT["axe"] == axe].set_index("n")
    out = {}
    previous_sheath = 0.

    def ramp(start, duration):
        return float(np.clip((thermal_time - start) / duration, 0., 1.))

    for n in range(1, n_metamers + 1):
        appeared = first_leaf + (n - 1) * phyllochron
        if thermal_time < appeared:
            break
        Ll, Gl, El = (float(dims.loc[n, name]) for name in ("Ll", "Gl", "El"))
        blade = Ll * ramp(appeared, 1.6 * phyllochron)
        sheath = Gl * ramp(appeared + phyllochron, phyllochron)
        internode = El * ramp(appeared + 2 * phyllochron, 1.5 * phyllochron)
        out[n] = {
            "blade": dict(length=blade, visible_length=max(0., blade - previous_sheath),
                          senesced_length=Ll * ramp(appeared + 5 * phyllochron, 3 * phyllochron) * (blade > 0),
                          age=thermal_time - appeared),
            "sheath": dict(length=sheath, visible_length=max(0., sheath - previous_sheath), senesced_length=0.),
            "internode": dict(length=internode, visible_length=0.5 * internode, senesced_length=0.),
        }
        out[n]["blade"]["senesced_length"] = min(out[n]["blade"]["senesced_length"], blade)
        previous_sheath = sheath
    return out


def triangles(geometry):
    """A PlantGL geometry's triangles in the scene's frame, as a (T, 3, 3) array."""
    tesselator = pgl.Tesselator()
    geometry.apply(tesselator)
    mesh = tesselator.result
    points = np.array([[p.x, p.y, p.z] for p in mesh.pointList], dtype=float)
    indices = np.array([[i[0], i[1], i[2]] for i in mesh.indexList], dtype=int)
    return points[indices] if indices.size else np.zeros((0, 3, 3))


def record(g) -> Record:
    vertices = [v for v in g.vertices() if v != g.root]
    names = [name for name in g.property_names() if name not in ("geometry", "edge_type", "label", "color")]
    props = {name: {v: value for v, value in g.property(name).items() if v in set(vertices)} for name in names}
    geometry = g.property("geometry")
    return Record(vertices=vertices, label={v: g.label(v) for v in vertices}, scale={v: g.scale(v) for v in vertices},
                  complex={v: g.complex(v) for v in vertices}, parent={v: g.parent(v) for v in vertices},
                  edge_type={v: g.edge_type(v) for v in vertices}, props=props,
                  triangles={v: triangles(geometry[v]) for v in vertices if v in geometry})


def legacy_run(n_steps=40, thermal_time_step=10., leaves="echap", seed=1, record_every=1, record_states=True,
               phyllochron=PHYLLOCHRON, plant_azimuth=0., inclinations=None, n_sect=1, plants=None):
    """*n_steps* steps of *thermal_time_step* °Cd of the legacy's dynamic geometry on the given lengths, the plant at
    the origin, its azimuth *plant_azimuth* (°; the legacy draws a random one). *inclinations*: {(axis, rank):
    degrees} of internodes (the legacy's AdelDyn sets none). *n_sect*: the number of elements of a blade. *plants*: a
    stand of several plants ([Planting]; then *plant_azimuth* is not used), the same tables for all (AdelDyn's
    add_metamer reads its phytomer table's first rows whatever the plant)."""
    axes, phytoT = harness_tables(leaves)
    plants = plants or [Planting(azimuth=plant_azimuth)]
    adel = AdelDyn(seed=seed, leaves=leaf_database(leaves), nsect=n_sect, nplants=len(plants))
    adel.positions = [(float(p.x), float(p.y), 0.) for p in plants]
    adel.plant_azimuths = np.array([float(p.azimuth) for p in plants])
    stand_axes = pd.concat([axes.assign(plant=number) for number in range(1, len(plants) + 1)], ignore_index=True)
    run = LegacyRun(adel=adel)
    g = adel.build_stand(stand_axes)
    labels = g.property("label")
    for step in range(n_steps + 1):
        for number, plant in enumerate(plants, start=1):
            lengths = given_lengths(start(axes) + step * thermal_time_step - plant.delay, axes, phytoT,
                                    phyllochron=phyllochron)
            for _, axis in axes.iterrows():
                name = axis["axe"]
                if not lengths[name]:
                    continue
                if name != "MS":
                    add_axe(g, name, plant_number=number, axis_properties=dict(
                        HS_final=float(axis["HS_final"]) if axis["HS_final"] == axis["HS_final"] else float(axis["nf"]),
                        nff=int(axis["nf"]), hasEar=bool(int(axis["hasEar"])), azimuth=0))
                _, _, metamers = find_metamers(g, number, name)
                for _ in range(len(lengths[name]) - (len(metamers) - 1)):            # metamer0 aside
                    adel.add_metamer(g, phytoT, plant=number, axe=name)
                _, _, metamers = find_metamers(g, number, name)
                for metamer in metamers:
                    rank = int(labels[metamer][len("metamer"):])
                    for organ in g.components(metamer):
                        values = lengths[name].get(rank, {}).get(labels[organ])
                        for property_name, value in (values or {}).items():
                            g.property(property_name)[organ] = value
                        if labels[organ] == "internode" and (name, rank) in (inclinations or {}):
                            g.property("inclination")[organ] = inclinations[(name, rank)]
        g = adel.update_geometry(g, phyllochron={"MS": phyllochron, "T1": phyllochron})
        if record_states and (step % record_every == 0 or step == n_steps):
            run.steps.append(record(g))
            run.recorded_steps.append(step)
    run.g = g
    return run


def turtle_placement(triangles, plant, axis):
    """The legacy's triangles (cm) of a plant's (*plant*: its Planting) element on *axis*, without update_geometry's
    second placement (the plant's position again, its azimuth's degrees read as radians) and T1 offset: as its turtle
    placed them, ShootArchitecture's placement."""
    angle = -plant.azimuth                             # the legacy gives the degrees to an angle in radians
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.], [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    out = (triangles - [plant.x, plant.y, 0.]) @ rotation.T
    if axis.startswith("T1"):
        out[..., 2] -= T1_OFFSET
    return out


def legacy_canopy(record, plants):
    """(vertices (T, 3, 3) in m, green (T,)) of a Record's elements, placed as the legacy's turtle placed them
    (turtle_placement; *plants*: the run's [Planting]); green: the element is more green than senesced."""
    vertices, green = [], []
    is_green = record.props.get("is_green", {})
    for element, triangles in record.triangles.items():
        if record.scale[element] != 5 or not len(triangles):
            continue
        axis = record.complex[record.complex[record.complex[element]]]
        plant = record.complex[axis]
        number = int(record.label[plant][len("plant"):])
        vertices.append(turtle_placement(triangles, plants[number - 1], record.label[axis]) * 1e-2)
        green.append(np.full(len(triangles), bool(is_green.get(element, True))))
    if not vertices:
        return np.zeros((0, 3, 3)), np.zeros(0, dtype=bool)
    return np.concatenate(vertices), np.concatenate(green)
