# OpenAlea.Adel

**Authors** : C. Fournier, C. Pradal, B. Andrieu

**Contributors** :

-   M Abichou (wheat parameteristion, documentation, plantgen),
-   C Chambon (plantgen, documentation, povray, fit, io widgets)

**Institutes** : INRAe, CIRAD

**Status** : R and Python package (+ L-system)

**License** : Cecill-C

# About

[![Last version](https://anaconda.org/openalea3/openalea.adel/badges/version.svg)](https://anaconda.org/OpenAlea3/openalea.adel/files)
[![Documentation Status](https://readthedocs.org/projects/adel/badge/?version=latest)](https://adel.readthedocs.io/en/latest/?badge=latest)
[![Licence](https://anaconda.org/openalea3/openalea.adel/badges/license.svg)](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html)
[![Platform](https://anaconda.org/openalea3/openalea.adel/badges/platforms.svg)](https://anaconda.org/openalea3/openalea.adel)
[![Downloads](https://anaconda.org/openalea3/openalea.adel/badges/downloads.svg)](https://anaconda.org/openalea3/openalea.adel)

## Description

OpenAlea.Adel (Architectural model of DEvelopment based on L-systems)
allows to simulate the 3D architectural development of the shoot of
gramineaous plant: where each leaf, sheath and stem piece of a plant is, how long and wide it is, and its shape, as
the plant grows. Light models then use these 3D plants, for instance to compute the light each leaf intercepts.

## Content

The package hosts generic data structure and simulation tools for
gramineous plants(Fournier & Pradal, unpublished), the Adel-Maize
(Fournier & Andrieu, 1998), Adel-Wheat (Fournier et al. 2003) models,
together with the wheat parameterization model of Abichou et al. (2013)
and the plastic leaf model of Fournier & Pradal (2012)


## Installation

### Users

```bash
mamba env create -n adel -c openalea3 -c conda-forge openalea.adel 
``` 

The metafspm port below (`ShootArchitecture`, `openalea.adel.simulation`) and its tutorials depend on
`openalea.metafspm`, which the developers' environment brings; they need the developers' installation until metafspm's
`release2026` and encoding3d are published.

### Developers

```bash
git clone -b release2026 https://github.com/openalea/metafspm.git metafspm
git clone -b other_platforms_bindings https://github.com/openalea/encoding3d.git encoding3d
git clone -b to_metafspm https://github.com/openalea/adel.git adel
cd adel
mamba env create -n adel_dev -f ./conda/environment.yml
mamba activate adel_dev
pip install -e ../encoding3d --no-deps                # until encoding3d is published
pip install -e ../metafspm --no-deps                  # until metafspm's release2026 is published
```

The environment brings `openalea.metafspm` (ShootArchitecture's framework, from `openalea3`), in an older version
until `release2026` is published: the two lines above install metafspm and encoding3d from their clones over it.

Tests: `pytest test` (the comparisons with the legacy need PlantGL; plantgen's need R through rpy2).

## Getting started

Grow a wheat plant for 80 days at 15 °C with the default parameters (Adel's sample plant), and look at it:

```python
from openalea.adel.simulation import simulate

result = simulate(days=80)
print(result.summary)                         # per day: thermal time, Haun stage, tillers, height, leaf area
result.draw()                                 # the plant in 3D
result = simulate({"phyllochron": 100.}, days=80)     # leaves appearing every 100 °Cd
result = simulate(days=80, outputs="outputs")         # and its files: run.json, scenario.json, summaries, frames.fspm.zip
```

`simulate(..., compute=Compute(...))` (`from openalea.metafspm import Compute`) sets how the run uses the machine:
the threads of numba's kernels, of the linear solves' blocks and of a reduction's kinds (`threads="auto"` by
default: the processors, each part taking those that pay, the results the same at any number; BLAS stays at one),
the device (`"cpu"`, or `"cuda"` with CuPy, experimental: it should work on suitable hardware but has not been tested on a GPU yet), the memory policy and
the graph systems' sparse solver (`linear_solver="superlu"` by default: one backend for every system, or a mapping
giving each its own, by its qualified name `"Component._system"`, else a property it declares such as `"symmetric"`,
else `"default"`; CHOLMOD, PARDISO and MUMPS optional, checked when the Scene is built). The default
`Compute()` suits a laptop; metafspm's user guide (`docs/user.md`, *The run configuration*) and *Running on a GPU*
(`docs/gpu.md`) give the details. The component's numba kernels (`openalea.adel.kernels`) compile at
their first call: `python -m openalea.metafspm warm_up` compiles them once, with metafspm's and the other ports'.

### User tutorials

Step by step, for newcomers to Adel, to plant models or to Python modelling, in `example/tutorials/user/`:

1. [Your first wheat plant](example/tutorials/user/01_first_run.ipynb): run the model, read and plot its results;
2. [Change a parameter](example/tutorials/user/02_change_a_parameter.ipynb): the parameters, another plant, tillers;
3. [Combine parameters](example/tutorials/user/03_combine_parameters.ipynb): two parameters together;
4. [Sensitivity analysis](example/tutorials/user/04_sensitivity_analysis.ipynb): which parameters matter most;
5. [Experimental data](example/tutorials/user/05_experimental_data.ipynb): measured dimensions and temperature,
   observations, fitting;
6. [With other models](example/tutorials/user/06_with_other_models.ipynb): a canopy lit by a light model.

## For developers

### Developer tutorials

In `example/tutorials/developer/`, after the first user tutorial:

0. [Overview](example/tutorials/developer/00_overview.ipynb): the modules, the original Adel, the validation;
1. [ShootArchitecture in metafspm](example/tutorials/developer/01_shoot_architecture_in_metafspm.ipynb): the
   component formalism, Adel's tables, a canopy, the tillers, the coupling contract, the comparison and the speed, a
   plant shared with other models (entities, an anchor);
2. [Extend the component](example/tutorials/developer/02_extend_the_component.ipynb): a new variable and process.

### ShootArchitecture: Adel's shoots in metafspm

`openalea.adel.architecture.ShootArchitecture` ports Adel's shoot architecture (its `AdelDyn`, as Wheat-BRIDGES uses
it) to a [metafspm](https://github.com/openalea/metafspm) `StructuralComponent`. It is R-free and computes its geometry
without PlantGL:

- the plants' topology in metafspm's MPG: axes (main stem and tillers), phytomers, organs (internode, sheath, blade),
  elements;
- the organs' dimensions and leaf shapes, from Adel's tables (`axeT`, `phytoT`) and leaf databases (`LeafShapes`:
  Adel's samples, echap's shapes by age; `LeafShapes.cnwheat()`, CN-Wheat's wheat: echap's Mercia profiles and
  Soissons' midribs);
- the elements, derived from the organs' lengths (hidden and visible parts, green and senesced, rolled and flat), and
  their areas;
- their placement by Adel's turtle (compiled), and their triangles (`triangles()`, metafspm's `Triangles`) for light
  models and display.

Development and senescence models write the organs' lengths. The coupling contract is in the module's documentation.
A canopy of plants is one population, in one MPG:

```python
from openalea.adel.development import random_delays, shoot_scene
from openalea.metafspm import stand
from openalea.adel.leaf_shapes import LeafShapes

scene = shoot_scene(stand(20, rows=4, plant_spacing=0.04, jitter=0.005, seed=0), delays=random_delays(20), leaves=LeafShapes.echap(), n_sect=2,
                    thermal_time_step=10.)
for _ in range(100):
    scene.run()                                        # 10 °Cd of development, then the shoot
architecture = scene.populations[0].instance.architecture
triangles = architecture.triangles()                   # every triangle and its element's node index
green = triangles.attributes["green"]                  # each triangle's green share (the senescence from the top)
```

- **Tables** (`openalea.adel.tables`): the tables given (`axes_table`, `phytomer_table`), or a plant generated from a
  few per-plant numbers (`leaf_number`, `tiller_number`, `longest_blade_length`, `longest_blade_width`: Adel's
  reference dimension profiles, scaled), or by default Adel's sample plant.
- **Development** (`openalea.adel.development`): `ScenarioDevelopment` is a simple deterministic scenario, linear in
  thermal time, with parameters for its timings; its thermal time sums the degree-days of the scene's air temperature
  (the forcing `air_temperature`, above `base_temperature`), or `thermal_time_step` a step without it. Adel-Wheat's R
  kinetics are not ported. `ScenarioShoot` is the plant model of both components (coupled by their variables' names,
  without a translator), and `shoot_scene` builds a scene of such plants.
- **Tillers:** vertical by default, as AdelDyn's. With `tiller_inclination=1` they are inclined by static Adel's rules
  (Adel.R's getIncT and axe_inclination: incT ramping with the Haun stage, straightening at dredT).
- **Ears:** on an axis with an ear (the axes table's `hasEar`) whose phytomer table has its ear rows (ranks nf + 1 to
  nf + 3, Adel's sample plant's main stem), its peduncle, ear and awns are organs of their own, cylinders of their
  rows' lengths and diameters, placed and lit as the internodes. `ScenarioDevelopment` brings them out
  `ear_start` (1.6) phyllochrons after the flag leaf, the peduncle growing for `peduncle_duration` (3).
- **Geometry options** (ShootArchitecture's settings, as Adel's): `twist` (°, the blades' twist about their midrib,
  0), `stem_mesh` (`"slim"`, 8 triangles a cylinder, or `"classic"`, PlantGL's 6-slice cylinder, 24) and `face_up`
  (False: the blades in the turtle's frame; True: in a frame whose z is the vertical).
- **Plants and the soil:** each plant is sown at minus its sowing depth (metafspm's `Planting`, m), by default 0.9
  times its coleoptile's length; T1 starts at the coleoptile's top, as T2.
- **Display:** metafspm's matplotlib drawings of the triangles (`openalea.metafspm.display`), with the shoot's
  colours (`openalea.adel.display.state_colors`).
- **A plant other models share** (a crown model's seed and subcrown internode, a root model's roots): ShootArchitecture
  declares its `entities` (`architecture.SHOOT_ELEMENTS`: its elements on its axes, those with a row in the axes
  table), and works on metafspm's `EntityView` of its elements only; on an anchored plant, the elements its main
  stem's first element bears (the subcrown's elements are StemElements too, on the crown's axis). Its steps per axis and organ, and ScenarioDevelopment's, leave the axes without a row in the axes table
  (`axis_row` -1) and their organs as they are. The plant's scenario entry `shoot_anchor` (a selector for metafspm's
  `MPG.find_anchor`, e.g. the crown's top `{"label": "StemElement", "last": True}`) hangs the main stem's first element
  `'<'` on that vertex, built by an initiator before ShootArchitecture; the main stem's base is then at the plant's x,
  y and at the anchor's `z2` (`shoot_base_height`), else the plant's z. Without an anchor, nothing changes
  (`test/test_anchored.py`).

Examples: `example/canopy.py`, `example/generated_plants.py`, `example/tillers.py`, `example/light_coupling.py` (a
toy light model coupled through the triangles). The docs' port pages (`doc/port/`) give the declarations, the equations,
the original's algorithms with every departure from them, and the performance.

#### The legacy and the comparisons

The code that is not ported (AdelWheat and its R kinetics, the MTG and L-system modules, the visual programming nodes,
...) is frozen in `openalea.adel.legacy`; the `wralea` entry point points there. `plantgen` stays at the top level.

`openalea.adel.comparison` runs the legacy's `AdelDyn` and the port on the same scenario. `test/comparison` compares
their topology, elements, anchor points and every triangle, on one plant and on a stand, to 1e-10 m.
`example/comparison` draws both canopies side by side and profiles them: the port is 4 to 12 times faster for 1 to
16 plants. The comparisons take the port's known differences into account:

- AdelDyn raises the T1 tillers' elements by 0.0285 (in its unit: 2.85 cm in Wheat-BRIDGES, which runs it in m); the
  port does not;
- AdelDyn places each plant twice (its turtle, then `transform_geom`, the azimuth's degrees read as radians); the port
  places it once, by the turtle;
- with several elements per blade, the legacy computes their bounds in cm and rounds them; where a bound falls on a
  point of the width profile, the mean widths can differ (in under 20 % of the elements, by under 10 %).
