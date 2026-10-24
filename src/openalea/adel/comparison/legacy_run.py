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

The scenario (given_lengths): on the main stem, metamer n appears at (n - 1) phyllochrons (add_metamer); its blade
grows linearly to its final length (phytoT's Ll) in 1.6 phyllochrons, its sheath (Gl) from 1 phyllochron for 1, its
internode (El) from 2 phyllochrons for 1.5; the blade visible above the previous sheath, the sheath above it too, the
internode half visible; the blade senesces from 5 phyllochrons for 3; its age is the thermal time since it appeared.
Lengths in the legacy's unit (cm, phytoT's), thermal time in °Cd.
"""
from dataclasses import dataclass, field

import numpy as np

import openalea.plantgl.all as pgl
from openalea.adel.legacy import data_samples
from openalea.adel.legacy.adel_dynamic import AdelDyn
from openalea.adel.legacy.echap_leaf import echap_leaves

PHYLLOCHRON = 88.92                                    # °Cd, the legacy's default for the main stem


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


def shape_indices(phytoT, leaves, nff, top_leaves=4):
    """phytoT with the leaf shape indices of the database *leaves*: the samples' as they are; echap's by its rule
    (echap_leaf.geoLeaf): class 2 for the top leaves (nf - n <= top_leaves), else 1, one shape per class (Lseed 1)."""
    if leaves != "echap":
        return phytoT
    adapted = phytoT.copy()
    adapted["Lindex"] = np.where(nff - adapted["n"] <= top_leaves, 2, 1)
    adapted["Lseed"] = 1
    return adapted


def given_lengths(thermal_time, n_metamers, phytoT, axe="MS", phyllochron=PHYLLOCHRON):
    """{rank n: {organ: {property: value}}} at *thermal_time* for the metamers appeared (rank n from 1)."""
    dims = phytoT.loc[phytoT["axe"] == axe].set_index("n")
    out = {}
    previous_sheath = 0.

    def ramp(start, duration):
        return float(np.clip((thermal_time - start) / duration, 0., 1.))

    for n in range(1, n_metamers + 1):
        appeared = (n - 1) * phyllochron
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
               phyllochron=PHYLLOCHRON):
    """*n_steps* steps of *thermal_time_step* °Cd of the legacy's dynamic geometry on the given lengths."""
    axeT, phytoT = data_samples.axeTable(), data_samples.phytoT()
    adel = AdelDyn(seed=seed, leaves=leaf_database(leaves))
    run = LegacyRun(adel=adel)
    g = adel.build_stand(axeT)
    nff = int(axeT.loc[axeT["axe"] == "MS", "nf"].iloc[0])
    phytoT = shape_indices(phytoT, leaves, nff)
    labels = g.property("label")
    for step in range(n_steps + 1):
        thermal_time = step * thermal_time_step
        lengths = given_lengths(thermal_time, nff, phytoT, phyllochron=phyllochron)
        metamers = [v for v in g.vertices(scale=3) if labels[v].startswith("metamer") and labels[v] != "metamer0"]
        for _ in range(len(lengths) - len(metamers)):
            adel.add_metamer(g, phytoT)
        metamers = sorted((v for v in g.vertices(scale=3) if labels[v].startswith("metamer")
                           and labels[v] != "metamer0"), key=lambda v: int(labels[v][7:]))
        for rank, metamer in enumerate(metamers, start=1):
            for organ in g.components(metamer):
                values = lengths[rank].get(labels[organ])
                if values:
                    for name, value in values.items():
                        g.property(name)[organ] = value
        g = adel.update_geometry(g, phyllochron={"MS": phyllochron, "T1": phyllochron})
        if record_states and (step % record_every == 0 or step == n_steps):
            run.steps.append(record(g))
            run.recorded_steps.append(step)
    run.g = g
    return run
