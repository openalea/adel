===========
Performance
===========

The times recorded while porting, on the comparisons' scenario (Adel's sample plant, organ lengths given at each step
of 10 °Cd). They were measured on a laptop often loaded by other work: they are indicative, each given with the
conditions it was measured in. The legacy's behaviours that the port changes are in :doc:`legacy`.


The legacy and the port, on stands
==================================

``example/comparison/profile_population.py`` grows stands of 1, 4 and 16 plants (echap's leaf shapes, positions on
a grid, azimuths and delays spread) for 60 steps, every step's geometry computed, as a light model coupled at each step needs it. Only the
models' calls are timed: on the legacy's side ``add_metamer`` and ``update_geometry`` (the harness' writing of the
given lengths, Python loops, left out); on the port's, the scene's step (``ScenarioDevelopment`` and
``ShootArchitecture``) and ``triangles()``. Numba's compilations are done before the timing. Recorded on 25 October
2026, the machine loaded (load average about 4.7):

====== ========= =============== ===================== =========
plants legacy    port: the step  port: ``triangles()`` speed-up
====== ========= =============== ===================== =========
1      1.86 s    0.37 s          0.10 s                3.9×
4      5.94 s    0.64 s          0.17 s                7.3×
16     12.80 s   0.79 s          0.29 s                11.9×
====== ========= =============== ===================== =========

The port's lead grows with the number of plants: its steps are vectorised over the population, in one MPG, where the
legacy visits every vertex in Python. For one plant, the port's step (6 ms) is mostly metafspm's cost per step; it
grows little with the plants.

The developer tutorial (``example/tutorials/developer/01_shoot_architecture_in_metafspm.ipynb``, its section 6)
runs the same profile for 1 and 4 plants. When the tutorial was written, it gave speed-ups of 2.7× (1 plant) and 5.5× (4 plants); its executed
version records:

====== ========= ============================ =========
plants legacy    port (step and triangles)    speed-up
====== ========= ============================ =========
1      1.81 s    0.76 s                       2.39×
4      3.40 s    0.52 s                       6.53×
====== ========= ============================ =========

The port's 4-plant run being shorter than its 1-plant run in that table shows how much these figures depend on the
machine's load at the time.


The blade meshes
================

On the harness' plant (378 elements, 1284 triangles), ``triangles()`` takes 20 ms when every blade element is meshed
anew and 1.1 ms when the meshes come from the cache (``ShootArchitecture._blade_triangles``: each blade element's local
mesh kept by MPG vertex while what it is made from is the same; only the grown and new elements are meshed again,
every element placed again). On the same plant, the legacy's turtle and meshes (``mtg_interpreter``) take 144 ms and
its whole ``update_geometry`` 307 ms.

The legacy's harness alone, the first measurement, takes about 1 s for 120 steps of 10 °Cd (12 metamers on the main
stem).


Speed-ups found by profiling
============================

Each keeps the legacy's results (the comparisons hold before and after):

- **the axes table read once** (``ShootArchitecture._axes_lookup``, cached): its pandas lookups were three quarters
  of the structural step;
- **the blade meshes compiled** (``kernels``, numba: the turtle's ``interpret`` and ``stem_inclinations`` and the
  meshes, in ``openalea.adel.kernels``; ``kernels.warm_up()``, run by metafspm's ``python -m openalea.metafspm
  warm_up``, compiles them once at install);
- **the midribs' simplification in plain arithmetic** on 3-vectors (``leaf_shapes._squared_distance``, in PlantGL's
  order of operations): fitting echap's database went from 2 s to 0.15 s, about 13× faster;
- **the default leaf shapes fitted once** (``architecture.default_leaves``): without a ``leaves`` parameter, Adel's
  samples database was fitted again at every call (``dict.setdefault`` evaluates its default argument before looking
  up the key); fitted once, a standalone run is about 50× faster.


Where the time goes
===================

- For one plant, metafspm's cost per step dominates the port's step (above).
- Before the axes table was read once, its pandas lookups were three quarters of the structural step.
- The turtle (``kernels.interpret``, compiled) places every element at each call; the blade meshes are rebuilt only for
  the elements that changed, so ``triangles()`` costs little once the canopy stops growing (1.1 ms against 20 ms on one
  plant).
- The traversal order (``turtle.traversal_order``) is computed once per change of topology by metafspm's
  ``data_structure.preorder`` (a compiled depth-first walk over the children's CSR, the coleoptile tillers deferred); only the
  children lists and the few deferred tillers are read in Python.


Sources
=======

The stand profile and the speed-ups found by profiling were recorded with the populations' comparison (Adel commit
90fb048); the blade mesh cache with the blade meshes (f771dc4); the default leaf shapes' fix with the standalone use
(c7792d4); the tutorial's figures in its executed notebook.
