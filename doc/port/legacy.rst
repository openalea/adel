=================================
The legacy: algorithms and quirks
=================================

ShootArchitecture ports Adel's dynamic shoot architecture, the R-free path through ``AdelDyn`` that Wheat-BRIDGES
uses (cnwgrass builds an ``AdelDyn``, calls ``add_metamer`` at each leaf emergence and ``update_geometry`` at every
step from the organ lengths its other models write, then reads the scene for Caribu). That path is frozen in
``openalea.adel.legacy``, verbatim except for its imports and data folders, with every module the port does not
replace (AdelWheat and its R kinetics, the MTG and L-system modules, the visual programming nodes, ...). This page
describes the legacy's algorithms as the port implements them, then registers every behaviour of the legacy that is
odd or that the port changes, with what the port does about it. The rules are in :doc:`equations`, the times in
:doc:`performance`.

The port does not reproduce the legacy bit for bit: PlantGL's and numpy's floating point operations differ. It is
compared with it within tolerances, on a deterministic scenario (``openalea.adel.comparison``), as listed at the end.


The legacy's chain and the port
===============================

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - step
     - legacy (``openalea.adel.legacy``)
     - port
   * - a stand of plants with an empty main stem
     - ``AdelDyn.build_stand``, ``mtg_editions.add_plant`` (a ``metamer0`` holding a collar)
     - ``ShootArchitecture.initiate_plant`` (the main stem's Axis, GrowthUnit, first phytomer); metafspm's planting
       table
   * - a new metamer, a new tiller
     - ``AdelDyn.add_metamer``, ``mtg_editions.add_vegetative_metamer``, ``add_axe`` (called by the development
       model)
     - ``ShootArchitecture._axis_phytomers`` (from the inputs ``axis_number``, ``phytomer_number``; metafspm's
       ``NewNodes``)
   * - the organs' dimensions
     - ``add_metamer`` reads ``phytoT``
     - ``tables.phytomer_dimensions``, ``ShootArchitecture._phytomer_values``; ``tables`` (given, generated, sample
       tables)
   * - the elements from the organs' lengths
     - ``mtg_editions.update_organ_elements``: ``newmtg.internode_elements``, ``sheath_elements``,
       ``blade_elements``, rebuilt at each step
     - ``ShootArchitecture._elements``, ``_blade_sectors`` (stable elements)
   * - the leaf shape database
     - ``geometric_elements.Leaves``, ``fitting.fit_leaves``, ``fit2``, ``simplify``, ``simplification.cost``;
       ``echap_leaf.echap_leaves``
     - ``leaf_shapes.LeafShapes`` (``samples``, ``echap``), ``fit_midrib``
   * - blade parts' areas and widths
     - ``Leaves.blade_elt_area``, ``Leaves.blade_elt_width``
     - ``LeafShapes.relative_area``, ``LeafShapes.mean_width``
   * - the turtle
     - ``mtg_interpreter.mtg_interpreter``, ``AdelVisitor``, ``AdelTurtle`` (PlantGL's ``PglTurtle``)
     - ``turtle.traversal_order``, ``kernels.interpret`` (compiled), metafspm's ``kernels`` frame operations
   * - the stem meshes
     - ``mtg_interpreter.compute_element``, ``StemElement_mesh``, ``slim_cylinder``
     - ``ShootArchitecture.triangles``, metafspm's ``slim_cylinders``, ``place``
   * - the blade meshes
     - ``Leaves.mesh``, ``incline_leaf``, ``fitting.mesh4``, ``leaf_element``, ``leaf_to_mesh``
     - ``kernels.incline``, ``blade_part``, ``strip_triangles``, ``blade_element_triangles``; a cache by vertex
   * - the plants' placement
     - the turtle from the plant's position, then ``update_geometry``'s ``transform_geom`` and the T1 offset
     - the turtle only (``_placement``)
   * - the scene for light models
     - ``Adel.scene`` (a PlantGL scene), triangulated by Caribu
     - ``ShootArchitecture.triangles()`` (metafspm's ``Triangles``: vertices and owner elements)
   * - the tillers' inclination
     - static Adel's R (``Adel.R``: ``getIncT``, ``axe_inclination``); none in ``AdelDyn``
     - ``turtle.base_inclination``, ``kernels.stem_inclinations``, with ``tiller_inclination=1``
   * - the organs' lengths
     - cnwgrass' morphogenesis and senescence (Wheat-BRIDGES); AdelWheat's R kinetics offline
     - any development model; ``development.ScenarioDevelopment`` (a linear scenario)


The structure
=============

**The MTG** (``mtg_editions``): five scales, plant / axis / metamer / organ / element. A plant starts with a main stem
holding a ``metamer0`` and its collar; each metamer holds an internode, a sheath and a blade, each organ a
``baseElement`` and a ``topElement`` sentinel between which its elements are inserted. A tiller (``add_axe``) is an
axis whose first metamer is branched on its bearer's metamer at its insertion rank, on its internode's base.
``add_metamer`` reads the metamer's dimensions from ``phytoT`` (the rows of the axis and the rank) and the shape key
from the leaf database (``get_leaf_key``).

**The port** (``ShootArchitecture``): the MPG's Plant / Axis / GrowthUnit / Phytomer / Organ / SubOrgan, one
GrowthUnit per axis; the elements created with their organ (an internode's and a sheath's ``HiddenElement`` and
``StemElement``, a blade's ``HiddenElement`` and ``n_sect`` ``LeafElement``\ s), never rebuilt; their lengths and
areas updated at each call. A legacy element is mapped to the port's by its organ and label; the legacy keeps only
the elements its organs' current lengths need, a subsequence of the port's.


The elements
============

``internode_elements`` and ``sheath_elements`` return a hidden element and a stem element (with ``split``, a green
and a senesced one); ``blade_elements`` a hidden element, then ``n_sect`` sectors along the mature shape, each with
its green, senesced, rolled and flat parts (the rolled part's diameter from the sector's area), its flat part's
position (srb, srt) and width. The areas of blade parts are Simpson's rule on the width profile's points within the
part (``blade_elt_area``), the widths the mean of the profile at the part's bounds and points strictly inside
(``blade_elt_width``). The port computes the same lengths, parts and widths, and the areas as exact integrals
(:doc:`equations`).


Leaf shapes
===========

``Leaves`` holds a database of midribs (``xydb``, by class and shape, by age class in a dynamic database) and width
profiles (``srdb``, by class), keyed by strings. ``fit_leaves`` fits each shape once (``fit3``: ``fit2``, then
``simplify``); a shape that fails is printed and discarded. ``fit2`` fits the midrib and the width profile with
scipy's smoothing splines; ``simplify`` simplifies the polyline of points (x, r, y) with PlantGL's ``Vector3`` by
``simplification.cost``, a heap of each interior point's squared distance to the line through its neighbours, and
rescales the width to keep the smooth shape's area. The port fits the same shapes with scipy and plain arithmetic on
3-vectors in PlantGL's order of operations, without PlantGL. ``get_age_index`` gives a dynamic database's age class
from the blade's age in phyllochrons minus 0.3.


The turtle and the meshes
=========================

``mtg_interpreter`` visits the MTG with openalea.mtg's ``pre_order2_with_filter`` (per plant, depth first, branches
before successors) and runs ``AdelVisitor`` at each vertex on an ``AdelTurtle``. The visitor keeps a context: the top,
the main stem's top and the tillers' base. At a plant's base it moves to the plant's position, sets the heading up and
the up towards −x, rolls by the plant's azimuth. A change of axis to a tiller's first metamer (``metamer0``, edge
``+``) saves the main stem's top after a main stem element, or returns to the tillers' base after another axis'; any
other change returns to the main stem's top (AD-37). A stem element's inclination turns the turtle down (an axis' first one, on a tiller,
rolled first by the axis' azimuth), the next ones after ``rollToVert`` and rolled back into their azimuth plane; its
azimuth rolls it. ``compute_element`` meshes a stem element of length at least ``min_length`` as a slim cylinder, a
blade element's rolled part as a slim cylinder and its flat part by ``Leaves.mesh`` (``incline_leaf``, x flipped,
``mesh4``: ``leaf_element``, ``leaf_to_mesh``, the triangles above ``min_length² / 100``), the flat part lifted by the
rolled length. ``AdelTurtle.transform`` places each mesh in the turtle's frame. The port's ``kernels.interpret`` follows
the same rules on arrays (compiled), and ``kernels`` the same mesh functions (compiled), in the legacy's order of
operations.


The quirk registry
==================

Each behaviour of the legacy that is odd, or that the port changes, with what the port does: **reproduced** (the
port computes what the legacy computes), **corrected** (the port computes the rule's exact result; named where an
option selects it), **fixed** (a bug, fixed in the port), **kept** (a modelling choice of the legacy kept), **changed**
(the port's design differs; the comparisons map one onto the other), **dropped** (not ported: files, options, dead or
out-of-scope code).

Structure and tables
--------------------

.. list-table::
   :header-rows: 1
   :widths: 8 52 40

   * - ID
     - legacy
     - port
   * - AD-01
     - ``find_label`` scans the graph on every call (quadratic in the number of vertices); properties carried by
       string labels
     - dropped (the MPG's indices; the steps vectorised over the population)
   * - AD-02
     - the elements removed and inserted again at each ``update_organ_elements``, by label (``StemElement``,
       ``LeafElement1`` … ``LeafElementn``)
     - changed: stable elements created with their organ; an element of zero length exists, with no area and no
       triangles
   * - AD-03
     - ``metamer0`` and its collar; the ``baseElement`` and ``topElement`` sentinels
     - dropped (the first phytomer is the first metamer; the comparisons leave the sentinels aside)
   * - AD-04
     - ``split=True``: separate green and senesced elements (``StemElementg``, ``LeafElement1s``, ...)
     - dropped: each element carries ``green_length``, ``green_area``, ``senesced_length``, ``senesced_area``
   * - AD-05
     - a blade has a ``HiddenElement`` (its part in the sheaths' tube), its area on the mature shape's base
     - reproduced
   * - AD-06
     - a tiller's first metamer branched on its bearer's metamer at its internode's base (``add_axe``); the static
       factory (``mtg_factory``) puts every tiller on the main stem
     - reproduced (``add_axe``'s insertion, at any bearer)
   * - AD-07
     - ``AdelDyn.add_metamer`` selects the ``phytoT`` rows by axis and rank only: the first matching row, whatever the
       plant
     - fixed: each plant reads its own rows (``table_plant``, the population's tables); the comparisons give every
       plant the tables' plant 1
   * - AD-08
     - the same property names at the organ and element scales (``length``, ``senesced_length``, ``area``)
     - changed: the organs' ``organ_length``, ``organ_senesced_length``, ``organ_area``; the elements' ``length``,
       ``senesced_length``, ``area``
   * - AD-09
     - lengths in cm (``conv_units``, ``scene_unit``), the unit conversions of ``AdelWheatDyn`` with ``.iteritems()``
       (Python 3 errors)
     - changed: m everywhere, the tables converted by ``length_unit`` (0.01) at their reading
   * - AD-10
     - ``new_mtg_factory`` unpacks three organs where an ear row gives one (a crash once ears appear)
     - dropped (``new_mtg_factory`` is not ported)
   * - AD-11
     - ears, peduncles and awns (geometry only)
     - dropped (not ported; the axes' ``has_ear`` kept as a descriptor)
   * - AD-12
     - the sample ``phytoT.csv`` has the main stem's ranks only
     - kept: the sample tables give each tiller the main stem's first ``nf`` ranks, lengths and diameters 0.9 times
       (``tables.with_tillers``), as the comparisons' harness does
   * - AD-13
     - plants defined only by tables (AdelWheat's ``setAdel`` or plantgen)
     - changed: a plant may be generated from four numbers by Adel's reference profiles (``tables.wheat_tables``), or
       be Adel's sample plant by default
   * - AD-14
     - the stand's positions drawn by an unseeded ``random``, the plants' azimuths by numpy's global seed
     - dropped: positions and rotations from metafspm's planting table; the comparisons give the legacy theirs
   * - AD-15
     - global state: the R environment, numpy's global seed
     - dropped (the port is R-free)
   * - AD-16
     - the development model calls ``add_metamer`` and writes the new metamer's lengths at once
     - changed: the organs a call creates get their lengths at the next call; ``ScenarioShoot`` runs its development
       and its shoot together (``CompositeModel.run_together``): the development's postsegmentation gives them after
       the shoot's segmentation
   * - AD-17
     - cnwgrass writes its own element values (the lengths of ``HiddenElement``, ``StemElement``, ``LeafElement1``)
       over Adel's after ``update_geometry``
     - changed: the elements are ShootArchitecture's; development models write the organs' lengths only
   * - AD-18
     - cnwgrass' morphogenesis writes the blades' ``shape_mature_length`` and ``shape_max_width``
     - reproduced: both are ShootArchitecture's, the phytomer table's ``Ll`` and ``Lw`` until a model gives a
       blade's ``final_blade_length`` and ``final_blade_width`` (its inputs)
   * - AD-53
     - the harness and cnwgrass add each axis by its name when its first leaf appears (``add_axe``: ``T1.1`` on
       T1's metamer 1, its bearers added first)
     - changed: a plant's ``axis_number`` counts its axes emerged, created in its axes table's order: the table must
       list them in their order of emergence (``emf1``), as Adel's tables do (``test_tiller_of_tiller`` sorts it)
   * - AD-54
     - a development adds a metamer per ``add_metamer`` call, one at a time
     - changed: a call creates every phytomer and axis its development asks for, each after the previous one
       (``ScenarioShoot.run_shoot`` segments again for the axes it created): one large step reaches the structure
       and the values of the small steps reaching the same thermal time, its elements created in another order

Elements and leaf shapes
------------------------

.. list-table::
   :header-rows: 1
   :widths: 8 52 40

   * - ID
     - legacy
     - port
   * - AD-19
     - the areas of blade parts by Simpson's rule on the width profile's points within each part (``blade_elt_area``),
       least exact on small parts (a tip just out)
     - corrected: the exact integral of the profile, linear between its points as every Adel function interpolates
       it. The legacy's areas differ by a median of 5e-4 to 7e-4, up to 7 % (samples) and 3.3 % (echap's)
   * - AD-20
     - widths: with several elements per blade, the sectors' bounds are computed in cm; where they fall on a point of
       the width profile, some round one ulp below (0.7499999999999999), and the point counts inside the element
       above
     - corrected: a point within 1e-9 of a bound counts as the bound (``BOUND_TOLERANCE``). With echap's profiles and
       4 elements per blade, 77 of 455 widths differ, by up to 9.6 %; the areas are not affected
   * - AD-21
     - ``leaf_keys`` sorts the classes' keys as strings: '10' before '2' (the wrong shapes from ten classes on)
     - corrected: classes in numerical order (the same for the databases here, of 6 and 2 classes)
   * - AD-22
     - the thresholds ``lvis > 1e-6`` (sectors computed) and ``lrolled < 1e-6`` (no rolled part) are absolute, in the
       legacy's unit (cm)
     - kept as numbers, now in m (1e-6 m)
   * - AD-23
     - a stem element is green when its green length is at least its senesced one (``>=``), a blade element when it
       is larger (``>``)
     - reproduced
   * - AD-24
     - the rolled part's diameter from the sector's area: :math:`S\,\rho/v / (\pi \rho)`, the area spread over the
       visible length as a cylinder's surface
     - reproduced
   * - AD-25
     - a shape that cannot be fitted is discarded with a printed message; the next shapes of its class take its
       number (``Lseed`` shifts); a dynamic shape is discarded when one of its age classes fails
     - reproduced (silently)
   * - AD-26
     - ``simplification.cost`` removes the second and the last but one points after its heap (the code's own comment:
       "bug in the algorithm")
     - reproduced
   * - AD-27
     - the simplification's distances in the space (x, r, y): the width mixed with the coordinates
     - reproduced
   * - AD-28
     - ``fit2``'s fallbacks when a spline fails (scipy's defaults, then order 1); ``simplify``'s width rescaled by the
       smooth shape's Simpson area over the polygon's trapezoids
     - reproduced (the fitted midribs equal the legacy's to 1e-12)
   * - AD-29
     - echap's shape classes (``geoLeaf``): class 2 where :math:`n_f - n \le 4`, else 1
     - reproduced (``tables.echap_classes``), for the sample and generated plants with a dynamic database
   * - AD-30
     - the age class of :math:`\text{age}/\varphi - 0.3`, :math:`\varphi` the main stem's phyllochron, and T1's for
       every tiller (``phyllochron.get(axe, phyllochron["T1"])``)
     - reproduced (``phyllochron``, ``tiller_phyllochron``: 88.92 and 91.6 °Cd, the legacy's defaults)
   * - AD-52
     - a blade without a shape (``Lindex`` or ``Lseed`` −999: its shape key None) keeps its sectors' lengths and flat
       parts, with no area, no width and no mesh (``blade_elt_area`` and ``compute_element`` without a key)
     - reproduced (``shape_class`` or ``shape_seed`` −1; NaN is read as −999 too), compared in
       ``test_tiller_of_tiller``

Placement
---------

.. list-table::
   :header-rows: 1
   :widths: 8 52 40

   * - ID
     - legacy
     - port
   * - AD-31
     - ``update_geometry`` raises the elements of the axes named T1… by 0.0285 in the scene's unit (0.0285 cm by
       default, 2.85 cm in Wheat-BRIDGES, which runs AdelDyn in m), a manual fix of the first tiller's position
     - dropped: T1 starts at the coleoptile's top (AD-33); the comparisons take the offset off the legacy's triangles
   * - AD-32
     - ``update_geometry`` places each plant a second time after the turtle (``transform_geom``): it translates the
       plant by its position again and rotates it by its azimuth, given in degrees, read as radians
     - fixed: each plant placed once, by the turtle. The comparisons undo the second placement
   * - AD-33
     - T1 (on the main stem's phytomer 1) starts at the plant's base, inside the coleoptile (phytomer 1's sheath)
     - corrected: the traversal visits T1 after the coleoptile's sheath and blade, so it starts at the coleoptile's
       top, where T2 starts, with the frame of its branching point; the MPG keeps it on phytomer 1. The comparisons
       move the legacy's T1 there (``port_run.coleoptile_top``)
   * - AD-34
     - plants at :math:`z = 0`
     - changed: a plant at minus its sowing depth (metafspm's planting convention), by default 0.9 times its
       coleoptile's final length (``shoot_scene(sowing_depth=)``); the comparisons sow at 0
   * - AD-35
     - ``AdelDyn`` sets every inclination to 0: the tillers grow vertical, inside the main stem's cylinder
     - kept by default; corrected with ``tiller_inclination=1``: static Adel's ``getIncT`` (incT ramping with the
       axis' Haun stage) and ``axe_inclination`` (the straightening at dredT), the azimuths azTb and azT
   * - AD-36
     - ``axe_inclination`` straightens a tiller where its stem reaches dredT horizontally; in the turtle the blades'
       rolled parts stack on the stem, which the rule does not count
     - kept: the straightened tillers reach dredT within 10 %, the excess up to 4 %
   * - AD-37
     - the turtle does not follow the topology: each element starts on the context's top; a tiller's first metamer
       (``metamer0``, the only one branched) after another axis restarts from the tillers' base, which an inclined
       tiller's first stem element saves after its azimuth roll; any other change of axis returns to the main stem's
       top saved. So a tiller of a tiller (T1.1) starts on the tillers' base, not on its bearer's phytomer, and its
       bearer continues after it on the main stem's top, not on its own, with its first inclined stem element (when
       still to come) treated as a next one: it loses its azimuth
     - reproduced by default (``kernels.interpret``: a tiller start only at the axis' first visit; T1's continuation
       on the main stem's top saved before its deferred start, AD-33), compared with T1.1 on T1's phytomer 1
       (``test_tiller_of_tiller``: anchors 1e-12 m, triangles 1e-10 m). Corrected with ``own_axis_top=1``: an axis
       continues on its own top where the turtle left it (a tiller with its first inclined element still to come),
       a tiller starts on the top its bearer's element left; T1 is then placed as without T1.1, and a plant without
       tillers of tillers as by Adel's turtle (tested with the given and static Adel's inclinations)
   * - AD-38
     - ``rollToVert`` with a vertical heading: PlantGL's convention, the left (0, −1, 0)
     - reproduced (metafspm's ``roll_to_vertical``, tested against ``PglTurtle``)
   * - AD-39
     - the turtle and every mesh recomputed over the whole canopy at each step
     - changed: every element placed at each call (compiled); the blade meshes kept by vertex while what they are made
       from is the same; the traversal recomputed only when the topology changes
   * - AD-40
     - the order of the traversal from openalea.mtg's ``pre_order2_with_filter``
     - reproduced (``turtle.traversal_order``), apart from T1 (AD-33)

Meshes
------

.. list-table::
   :header-rows: 1
   :widths: 8 52 40

   * - ID
     - legacy
     - port
   * - AD-41
     - a static database's blades (Adel's samples) are turned by ``incline_leaf`` with the blade's inclination as a
       factor of the shape's base angle; ``AdelDyn`` sets it to 0, so each shape's base is turned to the vertical.
       ``phytoT``'s blade angles (``incB``, ``dincB``, ``pAngle``, ``dpAngle``) are not read
     - kept (the blade's ``inclination``, 0; the angle columns not read)
   * - AD-42
     - a dynamic database's shapes carry their inclination (the factor 1)
     - reproduced
   * - AD-43
     - the shape's x flipped (``flipx``) to lie along the tiller turned by ``down``
     - reproduced
   * - AD-44
     - ``leaf_to_mesh``'s twist (``Leaves.twist``, 0 by default); ``classic`` stems (PlantGL's tapered cylinders, 6
       slices); ``face_up`` blades
     - dropped (no twist, slim cylinders, the turtle's frame)
   * - AD-45
     - one triangle at the tip when its width is below 0.001, an absolute value in the scene's unit
     - reproduced as 1e-5 m (Adel's 0.001 cm, ``kernels.TIP_WIDTH``)
   * - AD-46
     - the triangles of area at most ``min_length² / 100`` left out, ``min_length`` 0.01 cm
     - reproduced (``min_length`` 1e-4 m)
   * - AD-47
     - a blade element's rolled cylinder first, its flat part's mesh lifted by the rolled length
     - reproduced
   * - AD-55
     - a blade element's ``anchor_point`` is the turtle's frame, its blade's collar moved up by the rolled parts
       below it, for every element of the blade
     - kept as the frame's origin (``frame_x``, ``frame_y``, ``frame_z``, where it is drawn); its base and end
       (``x1`` … ``z2``, which other models read) along its blade: its rolled part's base or its flat part's midrib's
       first point, its flat part's midrib's last point
   * - AD-48
     - the scene as a PlantGL scene, triangulated by its user (Caribu's ``scene_to_cscene``), plants renumbered,
       translated and rotated point by point in the light model
     - changed: ``triangles()`` gives every triangle in the scene's frame with its element's node index; PlantGL is
       not used

Out of the ported path
----------------------

.. list-table::
   :header-rows: 1
   :widths: 8 52 40

   * - ID
     - legacy
     - port
   * - AD-49
     - Adel.R's ``visibility`` reads ``kin$Gl[i]`` and ``kin$Llvis[i]`` where ``w$Gt[i]`` is meant; ``setAdeluser``
       drops ``ssipars``; a partial ``$`` match (``ssisenT$rate`` for ``rate1``); ``openapprox`` defined twice with
       different behaviours
     - dropped (AdelWheat's kinetics are not ported)
   * - AD-50
     - ``run_adel`` uses ``nrem`` for ``nquot``; ``duplicated`` reads the property ``"position "``;
       ``plantgen_extensions``' ``tdisqp`` and ``decide_child_cohorts`` arity
     - dropped (in the legacy, not on the ported path)
   * - AD-51
     - duplicated implementations: ``mtg_factory`` and ``new_mtg_factory``, ``AdelDyn`` and ``AdelWheatDyn``,
       ``genGeoLeaf`` four times
     - dropped (one implementation of each rule)


The comparisons
===============

``openalea.adel.comparison.legacy_run`` drives the legacy's ``AdelDyn`` as cnwgrass does, on a deterministic scenario
of organ lengths given at each step on every axis (Adel's sample plant, the main stem and four tillers);
``port_run`` runs the port on the same scenario (``ScenarioShoot``). ``test/comparison``:

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - test
     - compared
     - within
   * - ``test_structure``
     - every axis' phytomers and organs at each step, the element layouts, dimensions, shapes, the tillers' insertions
     - equal (the legacy's elements a subsequence of the port's)
   * - ``test_elements``
     - every element's lengths (green, senesced, rolled, hidden), flat part, width; areas
     - 1e-9 (relative); areas within 8 %, median within 2e-3 (AD-19); widths: fewer than a fifth differ, by less
       than 10 % (AD-20)
   * - ``test_leaf_shapes``
     - the fitted midribs, the age classes (both databases)
     - 1e-12; equal
   * - ``test_geometry``
     - every drawn element's anchor point and triangles: echap's and the samples' shapes, 1, 3 and 4 elements per
       blade, a plant azimuth of 37° with inclined internodes
     - anchors 1e-12 m; triangles 1e-10 m, at the plant azimuth 0 (AD-31 and AD-33 undone)
   * - ``test_population``
     - four plants (positions, azimuths, delays) in one MPG against the legacy's stand, plant by plant
     - anchors and lengths 1e-12 m, triangles 1e-10 m
   * - ``test_tiller_inclination``
     - every stem organ's inclination and every axis' azimuths with ``tiller_inclination=1``, against a line by line
       transcription of ``getIncT`` and ``axe_inclination`` (R is not installed in the environment)
     - 1e-9°; dredT reached within 10 % (AD-36)
   * - ``test_tiller_of_tiller``
     - a tiller of a tiller (T1.1) and two blades without a shape: the organs, the elements, the anchor points and
       triangles (vertical, inclined, a plant azimuth of 37°); the correction ``own_axis_top``
     - as ``test_structure``, ``test_elements`` and ``test_geometry`` (AD-37, AD-52)
   * - ``test_scenario_development``
     - ``ScenarioDevelopment`` against the harness' given lengths, step by step
     - 1e-10


Constants
=========

============================================ ============================================================
``min_length``                               0.01 cm (1e-4 m): the shortest element drawn
the blades' triangles                        area above ``min_length² / 100``
the tip triangle                             width below 0.001 cm (1e-5 m): one triangle
sectors, rolled parts                        ``lvis > 1e-6``, ``lrolled < 1e-6`` → 0
the age class                                ``age / phyllochron − 0.3``; phyllochrons 88.92 (MS), 91.6 (T1, every tiller) °Cd
fitting                                      midrib spline smoothing 0.01, order 3, 100 points; profile 0.0001, 200 points
discretisation                               9 points (samples), 7 (echap's)
the T1 offset (dropped)                      0.0285 in the scene's unit
tillers' inclination (static Adel)           3° until Haun stage 1, then 30° per stage up to incT
the samples' tillers                         the main stem's ranks, 0.9 times
generated plants                             relative rank ``n − (N − 2.4)``; tillers ``(k + 0.6)`` phyllochrons, ``N − 2 − k``
                                             leaves, incT 60°, dredT 4 cm, azT 75°
the sowing depth                             0.9 times the coleoptile's final length
============================================ ============================================================
