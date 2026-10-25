"""
Plants from a few numbers, without tables: ShootArchitecture's per-plant parameters leaf_number, tiller_number,
longest_blade_length and longest_blade_width generate each plant's tables (openalea.adel.tables.wheat_tables: Adel's
reference dimension profiles, scaled). Three plants of one population, each its own, grown by ScenarioDevelopment
with inclined tillers; their blades' lengths along the main stem.

    python example/generated_plants.py [n_steps]     (generated_plants.png in the working directory)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

from openalea.adel.development import Planting, shoot_scene
from openalea.adel.display import draw

PLANTS = [Planting(x=0.), Planting(x=25.), Planting(x=50.)]
OWN = [dict(leaf_number=9, tiller_number=1, longest_blade_length=0.12, longest_blade_width=0.010),
       dict(leaf_number=11, tiller_number=2, longest_blade_length=0.18, longest_blade_width=0.014),
       dict(leaf_number=13, tiller_number=4, longest_blade_length=0.25, longest_blade_width=0.018)]

if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    scene = shoot_scene(PLANTS, plant_parameters=OWN, phyllochron=100., n_sect=2, tiller_inclination=1)
    for _ in range(n_steps + 1):
        scene.run()
    ds = scene.populations[0].data_structure
    architecture = scene.populations[0].instance.architecture

    figure = plt.figure(figsize=(14, 5))
    axes = figure.add_subplot(1, 2, 1, projection="3d")
    plant_of = ds.map_to(np.arange(len(PLANTS), dtype=float), "Plant", "node", "broadcast").astype(int)
    draw(axes, architecture.triangles(), colors=colormaps["tab10"](plant_of)[:, :3])
    axes.set_title(f"three generated plants, {n_steps * 10} °Cd")
    axes = figure.add_subplot(1, 2, 2)
    table_axes, phytomers, numbers, _ = architecture.tables()
    for k, (number, own) in enumerate(zip(numbers, OWN)):
        main = phytomers.loc[(phytomers["plant"] == number) & (phytomers["axe"] == "MS")]
        axes.plot(main["n"], main["Ll"], "o-", color=colormaps["tab10"](k),
                  label=f"{own['leaf_number']} leaves, {own['tiller_number']} tiller"
                        f"{'s' if own['tiller_number'] > 1 else ''}")
    axes.set_xlabel("rank on the main stem")
    axes.set_ylabel("blade length (cm)")
    axes.set_title("Adel's reference profile, scaled")
    axes.legend()
    figure.savefig("generated_plants.png", dpi=110)
    print(f"generated_plants.png: {len(PLANTS)} plants, their tables {sorted(set(numbers.tolist()))}")
