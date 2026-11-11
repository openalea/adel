# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `pyproject.toml` file following [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- [OpenAlea svg](./doc/_static/openalea_web.svg) logo for documentation.

### Fixed

- Fixed some scipy function call as well as pandas functions.
- Some documentation warnings.

### Changed

- The documentation configuration file `conf.py` with [OpenAlea Guidelines](https://openalea.readthedocs.io/en/latest/development/guidelines.html).
- Reformatted the formatting of file with [ruff](https://docs.astral.sh/ruff/).

### Removed

- ChangeLog.txt to use this CHANGELOG.md file instead.
- setup.py in favor of pyproject.toml.
- debian folder (not used)