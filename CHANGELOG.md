# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- metafspm fix-B (a packed run's frames folder removed once its archive is checked): `simulate(outputs=)` packs the
  frames into `frames.fspm.zip` when the run finishes (`SceneRecorder(pack_at_finish=True)`: the archive checked,
  then `frames/` removed); `example/exchange/export_frames.py` packs both of its exports and validates the archives,
  `import_frames.py` reads `outputs/frames.fspm.zip`; user tutorial 1 lists `frames.fspm.zip` (re-run).
- metafspm QG-3 (devplan_generic_kernels): metafspm's process pool and `Compute.workers` are removed; the README and
  `simulate`'s docstring describe `Compute(threads=)` alone (a reduction's kinds in its threads).

### Added

- `ShootArchitecture` (`openalea.adel.architecture`): Adel's shoot architecture as a metafspm StructuralComponent,
  R-free and without PlantGL: topology, elements, Adel's turtle (compiled), triangles for light models, optional
  tiller inclination by static Adel's rules; compared with the legacy's AdelDyn (`test/comparison`).
- `openalea.adel.tables` (sample and generated plants), `openalea.adel.leaf_shapes` (leaf databases fitted without
  PlantGL), `openalea.adel.development` (a simple scenario development, scenes of plants), `openalea.adel.display`
  (the shoot's colours for metafspm's drawings).
- Examples (`example/canopy.py`, `generated_plants.py`, `tillers.py`, `light_coupling.py`, `example/comparison`) and
  the tutorial `example/tutorials/tutorial_shoot_architecture.ipynb` (now
  `example/tutorials/developer/01_shoot_architecture_in_metafspm.ipynb`).
- `openalea.adel.simulation`, the users' entry point (`simulate`, `parameters`, `tables_from_measurements`); user
  tutorials (`example/tutorials/user/`: a first plant, a parameter, combinations, sensitivity analysis, experimental
  data, other models) and developer tutorials (`example/tutorials/developer/`); the docs split into a user guide and
  a developer guide.
- `LeafShapes.cnwheat()`: CN-Wheat's leaf shapes as a named preset (echap's Mercia widths, Soissons' midribs by leaf
  class and age).
- Adel's shoot on another model's crown: its own entities (`shoot_elements`), its coarse steps leaving other models'
  axes and organs unchanged; `shoot_anchor` (a `find_anchor` selector): the main stem's first element `<` on it, its
  base at the anchor's top (`shoot_base_height`); the crown's top found by `find_anchor(last=True)` (no `crown_top`
  mark).
- Examples of the exchange (`example/exchange`): `export_frames.py` writes the scene as encoding3d frames (formerly
  fspm-TF: the format other FSPM platforms and ParaView read), with a filtered export for logging beside the full
  one, and `import_frames.py` reads them back with openalea.encoding3d and with metafspm.

- `pyproject.toml` file following [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- [OpenAlea svg](./doc/_static/openalea_web.svg) logo for documentation.

### Fixed

- Fixed some scipy function call as well as pandas functions.
- `ShootArchitecture`: a call creates every phytomer and axis its development asks for (it crashed when a call asked
  for more than one phytomer per axis); `ScenarioShoot.grow` (was `run_shoot`) segments again for the axes a call created, one
  large step reaching the structure of small steps.
- `ShootArchitecture`: a blade element's `x1`...`z2` are its base and end along its blade (they were its blade's
  collar for all its elements); its frame's origin is `frame_x`, `frame_y`, `frame_z`.
- Some documentation warnings.

### Changed

- metafspm GK-2 (devplan_generic_kernels): `turtle.traversal_order(ds)` is metafspm's `ds.preorder` (a compiled
  depth-first walk, branches first, the coleoptile tillers deferred after their phytomer's last element) instead of
  a per-vertex walk with dict lookups; the same order and frame sources (test/comparison/test_traversal.py, against
  the former walk). Its argument is the DataStructure.
- metafspm P1 (the MPG's properties as columns): `ShootArchitecture`'s initiator writes the axes', phytomers' and
  elements' values with `g.set_values`, reads the plants' tables with `g.values_at`, and the turtle the edge types
  with `g.column("edge_type")`; `g.property(name)` is metafspm's read-only mapping now. The comparison with the
  legacy keeps openalea.mtg's API on the legacy's MTG.
- metafspm O-12 (part b): `ScenarioShoot` declares its parts and its `grow` step (`run_shoot` gone), nested by the
  plant models using it (`parts = {"shoot": ScenarioShoot, ...}`, `initiators = (..., ScenarioShoot, ...)`); the
  elements' `x1 ... z2` are topological (each element's base at its parent's end, its chord the heading times its
  length, a blade element's drawn chord), the turtle's ends their display ends `display_x1 ... display_z2`;
  `sown_stand` is gone (metafspm's `stand(n, rows=4, plant_spacing=0.04, jitter=0.005, seed=0)`) and
  `shoot_scene(depth=)` replaces `sowing_depth=` (metafspm's `Planting(depth=)`).

- `simulate(every="1 d")`: the summary's period is a duration (metafspm's `coupling.durations`: `"5 d"`), no longer a
  number of days; the frames every `"10 d"`.
- metafspm O-11: the triangles carry their green share (`architecture.triangles().attributes["green"]`): each
  triangle's area within its element's first `green_length` from its base (the senescence going down from the top),
  the distance along the element a cylinder's local z or a blade vertex's midrib abscissa (a second output of
  `kernels.blade_element_triangles`); the light models carry it along and the gas exchange's soil cut weighs by it
  (CNW-Grass's `emerged_links`: the green area above the soil). Behaviour change: an element Adel does not draw
  (shorter than `min_length`) has no triangles, so it now exchanges no gas.
- metafspm O-5: the numba kernels are in one module, `openalea.adel.kernels` (the turtle's `interpret` and
  `stem_inclinations`, the blade meshes: `blade_mesh.py` is gone); its `warm_up` is registered for metafspm's warm_up
  (entry point `openalea.metafspm.warm_up`) and run before the conda recipe's tests; metafspm's frame operations
  (`down`, `left_of`, `roll_right`, `roll_to_vertical`) come from `openalea.metafspm.kernels`.
- metafspm O-10 (forcings over the step): `ScenarioDevelopment`'s `air_temperature` is declared as its mean over the
  step (`forcing_variable`'s `at=` gone): the degree-days of a step are its mean temperature above the base times its
  duration (were its temperature at its start).
- metafspm O-1: `simulate(compute=)` and `shoot_scene(compute=)` forward metafspm's `Compute` (the Scene's run
  configuration: its threads, its memory trim; None: `Compute()`) to the Scene.
- The modules not ported to metafspm are frozen in `openalea.adel.legacy` (their tests and examples in `test/legacy`
  and `example/legacy`); the `wralea` entry point points there.

- The documentation configuration file `conf.py` with [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- Reformatted the formatting of file with [ruff](https://docs.astral.sh/ruff/).

### Removed

- ChangeLog.txt to use this CHANGELOG.md file instead.
- setup.py in favor of pyproject.toml.
- debian folder (not used)