"""
Importing fspm-TF files: the frames export_frames.py wrote (run first when missing), read two ways.

    python import_frames.py [--show]          figure in ./outputs/green_area_by_plant.png

1. As any platform reads them, with openalea.fspm_tf alone: each element's plant is found by following the tables'
   complex links up the scales (element → organ → phytomer → growth unit → axis → plant), which gives each plant's
   green area (leaves and stems), day by day, in the unit its declaration gives.
2. Back into metafspm, with openalea.metafspm.io.fspm_tf.import_frame: the last frame rebuilds Adel's DataStructure
   (its MPG, its variables, bit for bit), whose green area is the same.
"""
import argparse
import os

import matplotlib
import numpy as np
import pandas as pd

from openalea import fspm_tf
from openalea.metafspm.io.fspm_tf import import_frame

from export_frames import OUTPUTS, export

INDEX = os.path.join(OUTPUTS, "frames", "scene.fspm.json")


def plant_of(data_structure, set_name: str) -> np.ndarray:
    """The plant id of each entity of *set_name*, following the complex links up to the Plant entity set."""
    sets = {entity_set.name: entity_set for entity_set in data_structure.entity_sets}
    complex_of = {}                                                    # {id: its complex's id}, every scale
    for entity_set in data_structure.entity_sets:
        if "complex" in entity_set.columns:
            complex_of.update(zip(entity_set.columns["id"].tolist(), entity_set.columns["complex"].tolist()))
    plants = set(sets["Plant"].columns["id"].tolist())

    def up(identifier):
        while identifier not in plants:
            identifier = complex_of[identifier]
        return identifier
    return np.array([up(identifier) for identifier in sets[set_name].columns["id"].tolist()])


def green_area_by_day(index: str = INDEX) -> pd.DataFrame:
    """Each plant's green area (m2), one row per frame, read with openalea.fspm_tf."""
    rows = []
    for frame in fspm_tf.read_sequence(index):
        shoots = frame.data_structures[0]
        elements = shoots.entity_set("SubOrgan")
        assert elements.variables["green_area"].unit == "m2", elements.variables["green_area"].unit
        plants = plant_of(shoots, "SubOrgan")
        row = {"day": frame.time / 86400.}
        for number, plant in enumerate(shoots.entity_set("Plant").columns["id"]):
            row[f"plant {number + 1}"] = float(elements.columns["green_area"][plants == plant].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Adel's fspm-TF frames read back")
    parser.add_argument("--show", action="store_true")
    show = parser.parse_args().show
    if not os.path.exists(INDEX):
        export()
    table = green_area_by_day()
    print("Read with openalea.fspm_tf (green area per plant, leaves and stems, cm2):")
    print((table.set_index("day") * 1e4).round(2).to_string())

    last = fspm_tf.read_index(INDEX)["frames"][-1]["path"]
    identifier, ds = next(iter(import_frame(last).items()))            # {id: DataStructure}, as metafspm wrote them
    total = float(np.sum(ds.get("green_area")))
    print(f"\nRebuilt by metafspm from {os.path.basename(last)}: '{identifier}', a {type(ds).__name__} of "
          f"{ds.n_nodes()} nodes, {len(ds.entity_ids('Plant'))} plants; green area {total * 1e4:.2f} cm2")
    assert np.isclose(total, table.drop(columns="day").iloc[-1].sum(), rtol=0, atol=1e-15)

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(figsize=(6, 4))
    (table.set_index("day") * 1e4).plot(ax=axes, marker="o")
    axes.set_ylabel("green area, leaves and stems (cm²)")
    axes.set_title("Adel's wheat plants, read back from their fspm-TF frames")
    figure.tight_layout()
    path = os.path.join(OUTPUTS, "green_area_by_plant.png")
    figure.savefig(path, dpi=120)
    print(f"figure: {path}")
    if show:
        plt.show()


if __name__ == "__main__":
    main()
