"""
A canopy of 20 wheat plants grown by ShootArchitecture: four rows sown 15 cm apart, the plants a few days apart in
their development (GivenDevelopment, the harness' scenario, writes their organs' lengths), the leaf shapes by age
class (echap's). The canopy drawn at three times, its elements green or senesced, and its green and senesced areas
over time.

    python example/canopy.py [n_steps]               (canopy.png in the working directory)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

from openalea.adel.comparison.port_run import port_scene, sown_stand
from openalea.adel.display import draw, state_colors

N_PLANTS = 20
THERMAL_TIME_STEP = 10.                                # °Cd per step


def grow(n_steps, shown):
    """The canopy's triangles and its elements' state at the steps *shown*, and its areas at every step."""
    scene = port_scene("echap", plants=sown_stand(N_PLANTS), n_sect=2, thermal_time_step=THERMAL_TIME_STEP)
    ds = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture
    snapshots, areas = {}, []
    for step in range(n_steps + 1):
        scene.run()
        areas.append((np.sum(ds.get("green_area")), np.sum(ds.get("senesced_area"))))
        if step in shown:
            snapshots[step] = (architecture.triangles(), np.asarray(ds.get("is_green")).copy())
    return snapshots, np.array(areas)


if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    shown = (n_steps // 3, 2 * n_steps // 3, n_steps)
    snapshots, areas = grow(n_steps, shown)
    figure = plt.figure(figsize=(18, 5))
    for k, step in enumerate(shown):
        triangles, is_green = snapshots[step]
        axes = figure.add_subplot(1, 4, k + 1, projection="3d")
        draw(axes, triangles, colors=state_colors(is_green))
        axes.set_title(f"{step * THERMAL_TIME_STEP:.0f} °Cd, {len(triangles)} triangles")
    axes = figure.add_subplot(1, 4, 4)
    thermal_time = np.arange(n_steps + 1) * THERMAL_TIME_STEP
    axes.plot(thermal_time, areas[:, 0] * 1e4, color="tab:green", label="green")
    axes.plot(thermal_time, areas[:, 1] * 1e4, color="tab:olive", label="senesced")
    axes.set_xlabel("thermal time since the first leaf (°Cd)")
    axes.set_ylabel("area (cm²)")
    axes.set_title(f"the canopy's area, {N_PLANTS} plants")
    axes.legend()
    figure.subplots_adjust(left=0.02, right=0.98, wspace=0.25)
    figure.savefig("canopy.png", dpi=120)
    print(f"canopy.png: {N_PLANTS} plants, {n_steps} steps; green area at the end {areas[-1, 0] * 1e4:.0f} cm²")
