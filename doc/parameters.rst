=============================
The parameters in plain words
=============================

The main parameters of the plant model (``ShootArchitecture`` and ``ScenarioDevelopment``), with their default values.
``openalea.adel.simulation.parameters()`` lists every parameter with its value in the model's unit and its
description, and :doc:`port/components` gives the complete reference.

.. list-table::
   :header-rows: 1
   :widths: 26 52 22

   * - parameter
     - what it controls
     - default
   * - ``leaf_number``
     - number of leaves of the main stem
     - (12)
   * - ``tiller_number``
     - number of tillers
     - (4)
   * - ``longest_blade_length``
     - length of the longest blade
     - (19 cm)
   * - ``longest_blade_width``
     - width of the longest blade
     - (1.6 cm)
   * - ``phyllochron``
     - thermal time between two successive leaves
     - 88.92 °Cd
   * - ``blade_duration``
     - time a blade takes to grow, in phyllochrons
     - 1.6
   * - ``sheath_duration``
     - time a sheath takes to grow, in phyllochrons
     - 1
   * - ``internode_start``, ``internode_duration``
     - when an internode starts to grow and for how long, in phyllochrons
     - 2, 1.5
   * - ``senescence_start``
     - when a blade starts to senesce, in phyllochrons after its appearance
     - 5
   * - ``senescence_duration``
     - how long its senescence takes, in phyllochrons
     - 3
   * - ``tiller_inclination``
     - 0: the tillers grow vertical; 1: they lean, as in the field (static Adel's rules)
     - 0

The first four describe another plant than the default one. Not given (their default), the plant is Adel's sample
plant, with measured dimensions; given any of them, the plant is generated from the four, the missing ones taking the
values in brackets. A plant can also be given by Adel's two tables (``axes_table``, ``phytomer_table``), for example
from a main stem's measured dimensions (``openalea.adel.simulation.tables_from_measurements``).

The development follows the thermal time: ``simulate`` adds each day's mean temperature above a base temperature
(``base_temperature``, 0 °C by default) to it, in degree-days (°Cd).

To change a parameter: ``simulate({"phyllochron": 100.}, days=80)``.
