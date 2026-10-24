Adel's topology and geometry, frozen: the modules of `src/openalea/adel/` that `AdelDyn` (`adel_dynamic.py`),
`echap_leaf.py`, `Stand.py`, `data_samples.py` and `mtg_editions.py` import, from branch `to_metafspm` at 1b2f140,
verbatim. Only two kinds of lines are rewritten:

- their imports of each other (`openalea.adel.X` → `openalea.adel.legacy.X`);
- the data folders (`data_samples.py`, `geometric_elements.py`, `echap_leaf.py`): the package's (`../data`,
  `../echap_leaf_data`), not copied.

Lazy imports of modules outside this path (R through `AdelR`, `povray`, `leaf.curvature`, `wheat.extract_wheat`) still
point to the package's, as they are not part of the ported geometry.

It is the reference the port (ShootArchitecture) is compared with (`test/comparison/`, `example/comparison/`). It is
not developed further.
