"""
A wheat plant's tillers: their emergence, their insertion and their inclination. Each plant grows a main stem and
four tillers (T1 to T4, the axes table's rows), each inserted at its bearer's phytomer (T1 on the main stem's first)
when its first leaf appears (the table's emf1).

AdelDyn draws every axis vertical: the tillers grow inside the main stem. With tiller_inclination=1, ShootArchitecture
inclines them by static Adel's rules, from the axes table:
- the tiller's insertion azimuth (azTb), and its first phytomer's (azT);
- its first internode's inclination (incT), ramping with its Haun stage: 3° until stage 1, then 30° per stage;
- its straightening: once its stem reaches dredT horizontally from its insertion, it bends back to the vertical.

    python example/tillers.py [n_steps]              (tillers.png in the working directory)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

from openalea.adel.development import shoot_scene
from openalea.adel.display import draw
from openalea.adel.leaf_shapes import LeafShapes

THERMAL_TIME_STEP = 10.


def grow(n_steps, shown, tiller_inclination):
    """One plant: its triangles and their axes at the steps *shown*; every step's first internodes' inclinations."""
    scene = shoot_scene(leaves=LeafShapes.echap(xy_model="Soissons_byleafclass"), n_sect=2,
                        thermal_time_step=THERMAL_TIME_STEP, tiller_inclination=tiller_inclination)
    ds = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture
    snapshots, inclinations = {}, []
    for step in range(n_steps + 1):
        scene.run()
        first_internodes = (architecture.organ_names() == "internode") & (np.asarray(ds.get("organ_rank")) == 1)
        axis_of = ds.map_to(np.arange(len(ds.entity_ids("Axis")), dtype=float), "Axis", "Organ", "broadcast")
        inclination = dict(zip(architecture.axis_names()[axis_of[first_internodes].astype(int)],
                               np.asarray(ds.get("inclination"))[first_internodes]))
        inclinations.append(inclination)
        if step in shown:
            snapshots[step] = (architecture.triangles(), np.asarray(ds.owner("Axis")).copy(), architecture.axis_names())
    return snapshots, inclinations


def axis_colors(node_axis):
    return colormaps["tab10"](np.asarray(node_axis) % 10)[:, :3]


if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    shown = (30, 60, n_steps)
    inclined, inclinations = grow(n_steps, shown, tiller_inclination=1)
    vertical, _ = grow(n_steps, (n_steps,), tiller_inclination=0)

    figure = plt.figure(figsize=(16, 9))
    for k, step in enumerate(shown):
        triangles, node_axis, names = inclined[step]
        axes = figure.add_subplot(2, 3, k + 1, projection="3d")
        draw(axes, triangles, colors=axis_colors(node_axis))
        axes.set_title(f"{step * THERMAL_TIME_STEP:.0f} °Cd: {', '.join(names)}")
    triangles, node_axis, names = vertical[n_steps]
    axes = figure.add_subplot(2, 3, 4, projection="3d")
    draw(axes, triangles, colors=axis_colors(node_axis))
    axes.set_title(f"AdelDyn's vertical axes, {n_steps * THERMAL_TIME_STEP:.0f} °Cd")
    triangles, node_axis, names = inclined[n_steps]
    axes = figure.add_subplot(2, 3, 5, projection="3d")
    draw(axes, triangles, colors=axis_colors(node_axis), view=(90., -90.))
    axes.set_title("seen from above: the tillers' azimuths")
    axes = figure.add_subplot(2, 3, 6)
    thermal_time = np.arange(n_steps + 1) * THERMAL_TIME_STEP
    for k, name in enumerate(names):
        values = [step.get(name, np.nan) for step in inclinations]
        axes.plot(thermal_time, values, color=colormaps["tab10"](k), label=name)
    axes.set_xlabel("thermal time since the first leaf (°Cd)")
    axes.set_ylabel("the first internode's inclination (°)")
    axes.set_title("the tillers' inclination ramps with their Haun stage")
    axes.legend()
    figure.subplots_adjust(left=0.03, right=0.97, wspace=0.2, hspace=0.25)
    figure.savefig("tillers.png", dpi=110)
    print(f"tillers.png: the axes {', '.join(names)}")
