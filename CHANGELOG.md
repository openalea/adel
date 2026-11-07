# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
- The modules not ported to metafspm are frozen in `openalea.adel.legacy` (their tests and examples in `test/legacy`
  and `example/legacy`); the `wralea` entry point points there.

- The documentation configuration file `conf.py` with [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- Reformatted the formatting of file with [ruff](https://docs.astral.sh/ruff/).

### Removed

- ChangeLog.txt to use this CHANGELOG.md file instead.
- setup.py in favor of pyproject.toml.
- debian folder (not used)