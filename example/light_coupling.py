"""
ShootArchitecture coupled with a light model through its triangles (ShootArchitecture.triangles(), metafspm's
coupling.geometry.Triangles): the contract a port of Caribu will implement. The light model reads every triangle and
its element's node index, computes a light per triangle and gives it back per element (Triangles.per_owner).

Here a toy one, ToyLight: Beer's law in vertical columns, each triangle receiving exp(-k · the area of the triangles
above it in its column, per unit ground area). It is a component of the plant model, after the shoot: the canopy's
20 plants are one population, one MPG, so it sees all of them.

    python example/light_coupling.py [n_steps]       (light_coupling.png in the working directory)
"""
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from openalea.adel.architecture import SCALES, doc
from openalea.adel.comparison.given_development import GivenShoot
from openalea.adel.comparison.port_run import port_scene, plant_of, sown_stand
from openalea.adel.display import draw, value_colors
from openalea.metafspm.coupling.component import FunctionalComponent, input_variable, parameter, state_variable
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
    areas = triangles.areas()
    total = triangles.per_owner(areas, n_elements)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, triangles.per_owner(fraction, n_elements, weights=areas) / total, 0.)


@dataclass
class ToyLight(FunctionalComponent):
    """The light fraction each element receives, from the shoot's triangles; the light it absorbs, by its green
    area."""
    green_area: float = input_variable(**doc('m2', "The element's green area"), by='model_shoot_architecture',
                                       initialize=0., scale=SCALES.SubOrgan)
    light_fraction: float = state_variable(**doc('adim', "The light fraction the element receives"),
                                           state_variable_type="intensive", initialize=0., scale=SCALES.SubOrgan)
    absorbed: float = state_variable(**doc('m2', "The element's green area times its light fraction"),
                                     state_variable_type="extensive", initialize=0., scale=SCALES.SubOrgan)
    extinction: float = parameter(**doc('adim', "The extinction coefficient"), by="ToyLight", default=0.5)
    column_width: float = parameter(**doc('m', "The width of the columns"), by="ToyLight", default=0.05)

    shoot = None                                       # the StructuralComponent giving the triangles

    @rate
    def _light_fraction(self, green_area, extinction, column_width) -> tuple[np.ndarray, str, np.ndarray]:
        fraction = beer_law(self.shoot.triangles(), green_area.size, float(np.asarray(extinction).flat[0]),
                            float(np.asarray(column_width).flat[0]))
        return fraction, "absorbed", fraction * green_area


class LitShoot(GivenShoot):
    """GivenShoot, then the toy light on its triangles."""

    def __init__(self, data_structure, time_step):
        super().__init__(data_structure, time_step)
        self.light = ToyLight(data_structure=data_structure, time_step=time_step)
        self.light.shoot = self.architecture
        self.components.append(self.light)

    def run(self):
        super().run()
        self.light()


if __name__ == "__main__":
    n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    plants = sown_stand(20)
    scene = port_scene("echap", plants=plants, n_sect=2, thermal_time_step=THERMAL_TIME_STEP, model=LitShoot)
    ds = scene.populations[0].data_structure
    for _ in range(n_steps + 1):
        scene.run()
    architecture = scene.populations[0].instance.architecture
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
    absorbed = np.bincount(plant_of(ds, "node"), weights=np.asarray(ds.get("absorbed")), minlength=len(plants))
    print(f"light_coupling.png: {len(plants)} plants; the canopy's mean light fraction (green area weighted) "
          f"{absorbed.sum() / np.sum(ds.get('green_area')):.2f}; per plant, the most and least absorbed "
          f"{absorbed.max() * 1e4:.1f} and {absorbed.min() * 1e4:.1f} cm²")
