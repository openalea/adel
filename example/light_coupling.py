"""
ShootArchitecture coupled with a light model through its triangles (the DataStructure's triangles(), ShootArchitecture's:
metafspm's coupling.geometry.Triangles), the contract Caribu's port implements. The light model reads every triangle
and its element's node index, computes a light per triangle and gives it back per element
(Triangles.mean_per_owner, the area-weighted mean).

Here a toy one, ToyLight: Beer's law in vertical columns, each triangle receiving exp(-k · the area of the triangles
above it in its column, per unit ground area). It declares the geometry it reads (geometry_input) and the green area
(an input, ShootArchitecture's variable of the same name: no translator), and is coupled by couple() in the plant
model, after the shoot: the canopy's 20 plants are one population, one MPG, so it sees all of them. A light model
over several populations is an environment model on the plants' union instead (metafspm's EnvironmentModel, as
Caribu's CaribuLightModel).

    python example/light_coupling.py [n_steps]       (light_coupling.png in the working directory)
"""
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from openalea.adel.architecture import SCALES
from openalea.adel.development import ScenarioShoot, random_delays, shoot_scene
from openalea.metafspm.coupling.composite_wrapper import CompositeModel
from openalea.metafspm.scene.population import stand
from openalea.adel.leaf_shapes import LeafShapes
from openalea.metafspm.coupling.component import (FunctionalComponent, geometry_input, input_variable, setting,
                                                   state_variable)
from openalea.metafspm.display import draw, value_colors
from openalea.metafspm.solve.decorator import rate

THERMAL_TIME_STEP = 10.


def beer_law(triangles, n_elements, extinction, column_width):
    """Per element: the mean light fraction its triangles receive (area-weighted; 0 for an element without any),
    each triangle exp(-extinction · the projected area above it in its column / the column's ground area)."""
    if len(triangles) == 0:
        return np.zeros(n_elements)
    vertices = triangles.vertices
    centre = vertices.mean(axis=1)
    side, other = vertices[:, 1, :2] - vertices[:, 0, :2], vertices[:, 2, :2] - vertices[:, 0, :2]
    projected = 0.5 * np.abs(side[:, 0] * other[:, 1] - side[:, 1] * other[:, 0])
    _, column = np.unique(np.floor(centre[:, :2] / column_width).astype(np.int64), axis=0, return_inverse=True)
    column = column.reshape(-1)
    order = np.lexsort((-centre[:, 2], column))                     # by column, from the top
    cumulative = np.cumsum(projected[order])
    column_start = np.flatnonzero(np.r_[True, np.diff(column[order]) != 0])
    before_column = np.repeat(cumulative[column_start] - projected[order][column_start],
                              np.diff(np.r_[column_start, order.size]))
    above = np.empty_like(projected)
    above[order] = cumulative - projected[order] - before_column
    fraction = np.exp(-extinction * above / column_width ** 2)
    return triangles.mean_per_owner(fraction, n_elements)


@dataclass
class ToyLight(FunctionalComponent):
    """The light fraction each element receives, from the plants' triangles; the light it absorbs, by its green
    area."""
    geometry: object = geometry_input(description="The plants' elements' triangles (ds.triangles())")
    green_area: float = input_variable(unit='m2', description="The element's green area", by='model_shoot_architecture',
                                       default=0., scale=SCALES.SubOrgan)
    light_fraction: float = state_variable(unit='adim', description="The light fraction the element receives",
                                           kind="intensive", default=0., scale=SCALES.SubOrgan)
    absorbed: float = state_variable(unit='m2', description="The element's green area times its light fraction",
                                     kind="extensive", default=0., scale=SCALES.SubOrgan)
    extinction: float = setting(unit='adim', description="The canopy's extinction coefficient", default=0.5)
    column_width: float = setting(unit='m', description="The width of the columns", default=0.05)

    @rate
    def _light(self, green_area, extinction, column_width) -> dict:
        fraction = beer_law(self.data_structure.triangles(), green_area.size, extinction, column_width)
        return {"light_fraction": fraction, "absorbed": fraction * green_area}


class LitShoot(CompositeModel):
    """ScenarioShoot (nested), then the toy light on its triangles."""
    initiators = (ScenarioShoot,)
    parts = {"shoot": ScenarioShoot, "light": ToyLight}


if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    plants = stand(20, rows=4, plant_spacing=0.04, jitter=0.005, seed=0)
    scene = shoot_scene(plants, delays=random_delays(len(plants)), model=LitShoot,
                        leaves=LeafShapes.echap(xy_model="Soissons_byleafclass"), n_sect=2,
                        thermal_time_step=THERMAL_TIME_STEP)
    ds = scene.populations[0].data_structure
    for _ in range(n_steps + 1):
        scene.run()
    architecture = scene.populations[0].instance.shoot.architecture
    fraction = np.asarray(ds.get("light_fraction"))

    figure = plt.figure(figsize=(13, 5))
    axes = figure.add_subplot(1, 2, 1, projection="3d")
    draw(axes, architecture.triangles(), colors=value_colors(fraction, limits=(0., 1.)), shade=False)
    axes.set_title(f"the light fraction each element receives, {n_steps * THERMAL_TIME_STEP:.0f} °Cd")
    axes = figure.add_subplot(1, 2, 2)
    height = np.asarray(ds.get("z1"))
    lit = np.asarray(ds.get("green_area")) > 0
    axes.scatter(fraction[lit], height[lit] * 100, s=6, c="tab:green")
    axes.set_xlabel("the light fraction the element receives")
    axes.set_ylabel("the element's base height (cm)")
    axes.set_title("the light down the canopy (green elements)")
    figure.savefig("light_coupling.png", dpi=120)
    absorbed = ds.map_to(np.asarray(ds.get("absorbed")), "node", "Plant", "sum")            # per plant
    print(f"light_coupling.png: {len(plants)} plants; the canopy's mean light fraction (green area weighted) "
          f"{absorbed.sum() / np.sum(ds.get('green_area')):.2f}; per plant, the most and least absorbed "
          f"{absorbed.max() * 1e4:.1f} and {absorbed.min() * 1e4:.1f} cm²")
