"""
The time of ShootArchitecture's steps and triangles against the legacy's (AdelDyn's add_metamer and update_geometry),
on stands of 1 to 16 plants grown on the harness' scenario (the organs' lengths given at each step), every step's
geometry computed (as a light model at each step needs it).

Only the models' calls are timed: on the legacy's side, the harness' writing of the given lengths (Python loops) is
left out; on the port's, the scene's step includes ScenarioDevelopment's (vectorised, a small part of it).

    python example/comparison/profile_population.py [n_steps]
"""
import sys
import time

import numpy as np

from openalea.adel.comparison.legacy_run import given_lengths, harness_tables, leaf_database, start
from openalea.adel.comparison.port_run import port_scene
from openalea.adel.development import Planting
from openalea.adel.legacy.adel_dynamic import AdelDyn
from openalea.adel.legacy.mtg_editions import add_axe, find_metamers
import pandas as pd

THERMAL_TIME_STEP = 10.


def stand(n_plants, spacing=12.):
    """n_plants on a square grid (cm), their azimuths and delays spread."""
    side = int(np.ceil(np.sqrt(n_plants)))
    return [Planting(x=(k % side) * spacing, y=(k // side) * spacing, azimuth=(97. * k) % 360., delay=(13. * k) % 60.)
            for k in range(n_plants)]


def legacy_time(plants, n_steps, leaves="echap"):
    axes, phytomers = harness_tables(leaves)
    adel = AdelDyn(seed=1, leaves=leaf_database(leaves), nplants=len(plants))
    adel.positions = [(p.x, p.y, 0.) for p in plants]
    adel.plant_azimuths = np.array([p.azimuth for p in plants])
    g = adel.build_stand(pd.concat([axes.assign(plant=k) for k in range(1, len(plants) + 1)], ignore_index=True))
    labels = g.property("label")
    timed = 0.
    for step in range(n_steps + 1):
        for number, plant in enumerate(plants, start=1):
            lengths = given_lengths(start(axes) + step * THERMAL_TIME_STEP - plant.delay, axes, phytomers)
            for _, axis in axes.iterrows():
                name = axis["axe"]
                if not lengths[name]:
                    continue
                if name != "MS":
                    add_axe(g, name, plant_number=number, axis_properties=dict(
                        HS_final=float(axis["nf"]), nff=int(axis["nf"]), hasEar=False, azimuth=0))
                _, _, metamers = find_metamers(g, number, name)
                begin = time.perf_counter()
                for _ in range(len(lengths[name]) - (len(metamers) - 1)):
                    adel.add_metamer(g, phytomers, plant=number, axe=name)
                timed += time.perf_counter() - begin
                _, _, metamers = find_metamers(g, number, name)
                for metamer in metamers:
                    rank = int(labels[metamer][len("metamer"):])
                    for organ in g.components(metamer):
                        for property_name, value in (lengths[name].get(rank, {}).get(labels[organ]) or {}).items():
                            g.property(property_name)[organ] = value
        begin = time.perf_counter()
        g = adel.update_geometry(g)
        timed += time.perf_counter() - begin
    return timed, sum(1 for v in g.property("geometry"))


def port_time(plants, n_steps, leaves="echap"):
    scene = port_scene(leaves, plants=plants, thermal_time_step=THERMAL_TIME_STEP)
    architecture = scene.populations[0].instance.architecture
    stepping, meshing = 0., 0.
    for step in range(n_steps + 1):
        begin = time.perf_counter()
        scene.run()
        stepping += time.perf_counter() - begin
        begin = time.perf_counter()
        triangles = architecture.triangles()
        meshing += time.perf_counter() - begin
    return stepping, meshing, len(np.unique(triangles.owner))


if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    port_time(stand(1), 2)                                          # compilations aside
    print(f"{n_steps} steps of {THERMAL_TIME_STEP} °Cd, every step's geometry")
    print(f"{'plants':>6} {'legacy (s)':>11} {'port step (s)':>14} {'port triangles (s)':>19} {'speed-up':>9}"
          f" {'drawn elements (legacy, port)':>30}")
    for n_plants in (1, 4, 16):
        plants = stand(n_plants)
        legacy, legacy_drawn = legacy_time(plants, n_steps)
        stepping, meshing, port_drawn = port_time(plants, n_steps)
        print(f"{n_plants:>6} {legacy:>11.2f} {stepping:>14.2f} {meshing:>19.2f} {legacy / (stepping + meshing):>8.1f}x"
              f" {legacy_drawn:>15} {port_drawn:>14}")
