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
  (matplotlib drawings).
- Examples (`example/canopy.py`, `generated_plants.py`, `tillers.py`, `light_coupling.py`, `example/comparison`) and
  the tutorial `example/notebook/tutorial_shoot_architecture.ipynb`.

- `pyproject.toml` file following [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- [OpenAlea svg](./doc/_static/openalea_web.svg) logo for documentation.

### Fixed

- Fixed some scipy function call as well as pandas functions.
- Some documentation warnings.

### Changed

- The modules not ported to metafspm are frozen in `openalea.adel.legacy` (their tests and examples in `test/legacy`
  and `example/legacy`); the `wralea` entry point points there.

- The documentation configuration file `conf.py` with [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- Reformatted the formatting of file with [ruff](https://docs.astral.sh/ruff/).

### Removed

- ChangeLog.txt to use this CHANGELOG.md file instead.
- setup.py in favor of pyproject.toml.
- debian folder (not used)