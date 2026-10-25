"""
The legacy's canopy (AdelDyn) and ShootArchitecture's, side by side: four plants grown on the harness' scenario (their
positions, azimuths and delays differ), drawn at the end. The legacy's triangles are drawn as its turtle placed them
(legacy_run.turtle_placement: update_geometry's second placement and T1 offset undone), where test_population finds
them equal to the port's.

    python example/comparison/side_by_side.py [n_steps]      (side_by_side.png in the working directory)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

from openalea.adel.comparison.legacy_run import Planting, legacy_canopy, legacy_run
from openalea.adel.comparison.port_run import port_run
from openalea.adel.display import draw, state_colors
from openalea.metafspm.coupling.geometry import Triangles

PLANTS = [Planting(0., 0., 0., 0.), Planting(15., 0., 40., 30.), Planting(0., 12., 130., 60.),
          Planting(15., 12., 250., 15.)]

if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    options = dict(n_steps=n_steps, leaves="echap", n_sect=2, plants=PLANTS, record_every=n_steps)
    legacy = legacy_run(**options)
    port = port_run(**options)
    legacy_vertices, legacy_green = legacy_canopy(legacy.steps[-1], PLANTS)
    ds = port.scene.populations[0].data_structure
    port_triangles = port.scene.populations[0].instance.architecture.triangles()

    figure = plt.figure(figsize=(13, 6))
    for k, (title, vertices, colors) in enumerate((
            ("the legacy (AdelDyn)", legacy_vertices, state_colors(legacy_green)),
            ("the port (ShootArchitecture)", port_triangles.vertices,
             state_colors(np.asarray(ds.get("is_green")))[port_triangles.owner]))):
        axes = figure.add_subplot(1, 2, k + 1, projection="3d")
        draw(axes, vertices, triangle_colors=colors)
        area = Triangles(vertices, np.zeros(len(vertices))).areas().sum()
        axes.set_title(f"{title}\n{len(vertices)} triangles, {area * 1e4:.1f} cm²")
        print(f"{title}: {len(vertices)} triangles, their area {area * 1e4:.6f} cm²")
    figure.savefig("side_by_side.png", dpi=120)
