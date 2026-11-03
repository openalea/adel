=====================
Rules and equations
=====================

The rules of ShootArchitecture, as the port computes them, with the functions that compute them. Adel states most of
them as code (its tables' columns, its elements' functions, its turtle's context): they are written here. The
legacy's own algorithms, and where the port departs from them, are in :doc:`legacy`.

Notation: lengths in m, angles in degrees, thermal time in °Cd. An organ has a length :math:`l`, a visible length
:math:`l_v` (out of the sheaths' tube), a senesced length :math:`l_s` and, a blade, a rolled length :math:`l_r`; a
stem organ (internode, sheath) has a diameter :math:`d`; a blade has a mature length :math:`L`
(``shape_mature_length``), a maximal width :math:`W` (``shape_max_width``) and a shape class whose width profile
:math:`r(s)` gives the width relative to :math:`W` at the relative abscissa :math:`s \in [0, 1]` from the blade's base
to its tip. An axis' final leaf number is :math:`n_f`, a phytomer's rank :math:`n` (from 1).


The organs' dimensions (``tables``, ``architecture``)
=====================================================

**The tables.** A plant's dimensions come from Adel's two tables: the axes table (axeT: one row per axis, columns
``plant``, ``axe``, ``nf``, ``HS_final``, ``hasEar``, ``emf1``, optionally ``incT``, ``dredT``, ``azTb``, ``azT``) and
the phytomer table (phytoT: one row per plant, axis and rank ``n``). A plant's tables are, in order of precedence
(``tables.plant_tables``): the tables given (``axes_table``, ``phytomer_table``, the plant's rows by its
``table_plant``); tables generated when any of the per-plant numbers ``leaf_number``, ``tiller_number``,
``longest_blade_length``, ``longest_blade_width`` is given (``tables.wheat_tables``); Adel's sample plant
(``tables.sample_tables``). Each plant's tables are its ``shoot_tables`` (a Plant state, written at its initiation,
which the development reads as an input); ``tables.population_tables`` gathers them into one pair of tables for the
population, each distinct pair of tables (a source) once, its plants renumbered :math:`1000\,k + \text{plant}` for the
:math:`k`-th source, each plant's axis rows at an offset (written at its initiation, from the plants built before it).

**A phytomer** (``tables.phytomer_dimensions``, ``ShootArchitecture._phytomer_values``), at its creation, from its axis' row and its
rank, with :math:`u` its plant's tables' length unit (``length_unit``, Adel's cm: 0.01):

==================================== =========================== ========================================
variable                             organ                       value
==================================== =========================== ========================================
``shape_mature_length``              blade                       :math:`u\,\text{Ll}`
``shape_max_width``                  blade                       :math:`u\,\text{Lw}`
``diameter``                         internode, sheath           :math:`u\,\text{Ed}`, :math:`u\,\text{Gd}`
``organ_azimuth``                    internode                   ``Azim``
``shape_class``, ``shape_seed``      blade                       ``Lindex``, ``Lseed`` (−1 when −999)
``ntop``                             internode, sheath, blade    :math:`n_f - n + 1` (1: the flag leaf)
==================================== =========================== ========================================

The sheath's and internode's lengths ``Gl`` and ``El`` are read by a development model, not by ShootArchitecture.
``shape_mature_length`` and ``shape_max_width`` are ShootArchitecture's: a development model giving a blade's final
dimensions writes its inputs ``final_blade_length`` and ``final_blade_width`` (NaN: not given), which they take from
then on (``ShootArchitecture._elements``).

**An axis** (``ShootArchitecture._axis_values``), from its row: ``final_leaf_number`` :math:`= n_f`;
``final_haun_stage`` = ``HS_final`` (:math:`n_f` when missing); ``has_ear`` = ``hasEar``; ``table_inclination`` =
``incT``, ``straightening_distance`` :math:`= u\,\text{dredT}`, ``insertion_azimuth`` = ``azTb``, ``first_azimuth``
= ``azT`` (0 when the column or the value is missing). Its bearer and insertion rank come from its name
(``architecture.bearer_of``): ``MS`` has none; ``Tk`` is borne on the main stem's phytomer :math:`k`; ``T1.2`` on
T1's phytomer 2.

**Adel's sample plant** (``tables.sample_tables``, ``tables.with_tillers``): ``data/axeTable.csv`` (a main stem and
four tillers) and ``data/phytoT.csv`` (the main stem's ranks). A tiller with no rows of its own takes its plant's main
stem's first :math:`n_f` ranks, with ``Ll``, ``Lw``, ``Gl``, ``Gd``, ``El``, ``Ed`` multiplied by 0.9.

**A generated plant** (``tables.wheat_tables``), from :math:`N` leaves (``leaf_number``, 12 when not given),
:math:`K` tillers (4), the longest blade's length :math:`L_\text{max}` (0.19 m) and width :math:`W_\text{max}`
(0.016 m), and the phyllochron :math:`\varphi`. The main stem's ranks follow Adel's reference profiles
(``data/Wheat_dimension_profiles.csv``, plantgen's WheatDimensions) along the relative rank

.. math::

   x_n = \operatorname{clip}\big(n - (N - 2.4),\ x_\text{min},\ x_\text{max}\big),

the profiles' columns interpolated linearly at :math:`x_n`, then scaled:

.. math::

   \text{Ll}, \text{Gl}, \text{El} = \lambda\,(L_\text{blade}, L_\text{sheath}, L_\text{internode})(x_n), \qquad
   \lambda = \frac{L_\text{max}/u}{\max_n L_\text{blade}(x_n)},

.. math::

   \text{Lw}, \text{Gd}, \text{Ed} = \omega\,(W_\text{blade}, W_\text{sheath}, W_\text{internode})(x_n), \qquad
   \omega = \frac{W_\text{max}/u}{\max_n W_\text{blade}(x_n)},

with ``Azim`` 180, ``Lseed`` 1 and ``Lindex`` :math:`\min(n_\text{top}, C)` (:math:`C` the database's number of
classes), or, for a dynamic database, echap's classes (``tables.echap_classes``: 2 where :math:`N - n \le 4`, else 1).
The main stem has ``emf1`` 0 and no inclination. Tiller :math:`T_k` (:math:`k = 1 \dots K`, while
:math:`N - 2 - k \ge 1`) has :math:`n_f = N - 2 - k` leaves, ``emf1`` :math:`(k + 0.6)\,\varphi`, ``incT`` 60°,
``dredT`` 4 cm, ``azT`` 75°, ``azTb`` 0°, and its ranks the main stem's times 0.9 (``with_tillers``). The tables are
cached: equal numbers give the same tables.

**The sowing depth** (``development.coleoptile_sowing_depth``): a plant of ``shoot_scene`` is placed at
:math:`z = -\delta`, by default :math:`\delta = 0.9\,u\,\text{Gl}(\text{MS}, 1)`, nine tenths of its coleoptile's
final length, so that the coleoptile's tip emerges.


The structure (``ShootArchitecture._axis_phytomers``)
=====================================================

At each call, an axis is given phytomers up to :math:`\lfloor \texttt{phytomer\_number} \rfloor`, and a plant axes up
to :math:`\lfloor \texttt{axis\_number} \rfloor`, in its axes table's order. A new axis comes with its first phytomer,
branched (``+``) on the first element of its bearer's internode at its insertion rank; while that phytomer does not
exist, the axis and the plant's axes after it wait. A phytomer holds an internode (a ``HiddenElement`` then a
``StemElement``), a sheath (the same, branched on the internode's last element) and a blade (a ``HiddenElement`` then
``n_sect`` ``LeafElement``\ s, on the sheath's last element) (``architecture.phytomer_layout``). All are created in
one batch per call; their elements are never rebuilt.

**On a plant other models share** (``architecture.shoot_elements``, the initiator's ``shoot_anchor``). With the plant's
scenario entry ``shoot_anchor``, a selector of metafspm's ``MPG.find_anchor`` (``label``, ``scale``, ``where``,
``first`` or ``last``), the main stem's first element follows (``<``) that vertex, built by an earlier initiator (a
crown's top StemElement, the last of them: ``{"label": "StemElement", "last": True}``), instead of starting the plant's
topology; it is the plant's ``shoot_base`` (an MPG property). The architecture's entities (its ``entities``:
metafspm's ``EntityView``) are the elements (``HiddenElement``, ``StemElement``, ``LeafElement``) and, on an anchored
plant, only those its ``shoot_base`` bears, itself included: the nodes whose path from the root meets it. The
selection uses what a new element has at its creation, its label and its parents. The coarse locations are shared: an
axis without a row in the axes table (``axis_row`` :math:`= -1`, another model's) and its organs keep their values in the steps per
axis and per organ (``_axis_phytomers``, ``_first_inclination``, ``_elements``' organ areas,
``ScenarioDevelopment``'s), and its name is ``""``. The traversal starts each plant at its elements without a parent
among them. The main stem's base is at :math:`(x, y, z_a)`, :math:`z_a` the anchor's ``z2`` when the model that built it
gave it (``shoot_base_height``, written at initiation), else the plant's :math:`z`: a plant sown at the depth of the
crown's top grows the same elements.


The elements (``ShootArchitecture._elements``)
==============================================

**Stem organs** (internodes and sheaths): a cylinder of the organ's diameter. With the hidden length
:math:`l_h = \max(l - l_v, 0)`, the green length :math:`l_g = l_v - \min(l_s, l_v)` and the visible senesced length
:math:`l_{vs} = l_v - l_g`:

.. math::

   \text{hidden element:}\quad \ell = l_h,\ A = \pi\, l_h\, d; \qquad
   \text{stem element:}\quad \ell = l_v,\ A = \pi\, l_v\, d,\ A_g = \pi\, l_g\, d,\ A_s = \pi\, l_{vs}\, d,

the stem element green when :math:`l_g \ge l_{vs}`.

**Blades.** The visible part of a growing blade is the top of its mature shape: the blade's abscissa :math:`\sigma`
(in m from the mature shape's base, :math:`0 \le \sigma \le L`) is visible above :math:`L - l_v`. The hidden element
has :math:`\ell = l_h` and, when :math:`l_h > 0`,

.. math::

   A_h = L\,W \int_{(L - l)/L}^{(L - l_v)/L} r(s)\,ds .

The rolled part is :math:`l_r' = \operatorname{clip}(l_r, 0, l_v)`, 0 below 1e-6, the flat part
:math:`l_f = l_v - l_r'`, the senesced part :math:`l_{vs} = \min(l_s, l_v)`. With :math:`n` sectors
(``n_sect``, ``_blade_sectors``) of length :math:`\Delta = L/n`, sector :math:`k` covers
:math:`[(k - 1)\Delta, k\Delta]` of the mature shape; with :math:`t_k = k\Delta`:

.. math::

   v_k = \min\big(\Delta,\ \max(0,\ t_k - (L - l_v))\big), \qquad b_k = t_k - v_k, \qquad
   g_k = \min\big(t_k,\ \max(b_k,\ L - l_{vs})\big),

its visible length :math:`v_k`, its green part :math:`[b_k, g_k]` and its senesced part :math:`[g_k, t_k]`:

.. math::

   A_{g,k} = L\,W \int_{b_k/L}^{g_k/L} r(s)\,ds, \qquad A_{s,k} = L\,W \int_{g_k/L}^{t_k/L} r(s)\,ds, \qquad
   A_k = A_{g,k} + A_{s,k},

its length :math:`v_k`, green when :math:`g_k - b_k > t_k - g_k`. Its flat part and rolled part:

.. math::

   f_k = \min\big(v_k,\ \max(0,\ t_k - (L - l_f))\big), \qquad \rho_k = v_k - f_k, \qquad
   d_{\rho,k} = \frac{A_k\,\rho_k / v_k}{\pi\,\rho_k} = \frac{A_k}{\pi\, v_k} \quad (\rho_k > 0),

the rolled part a cylinder of diameter :math:`d_{\rho,k}` (``rolled_diameter``), the flat part's position relative to
the visible length (Adel's srb and srt; NaN without a flat part):

.. math::

   \text{flat\_base} = \frac{t_k - f_k - (L - l_v)}{l_v}, \qquad \text{flat\_top} = \frac{t_k - (L - l_v)}{l_v}.

Its ``width`` is the larger of its green and senesced parts' mean widths (``LeafShapes.mean_width``, below), times
:math:`W`. Sectors are computed for a blade with :math:`L > 0` and :math:`l_v > 10^{-6}`; without a shape class, with no area
and no width (Adel's, its shape key None), and no mesh.

**The organs' visible areas** (``organ_area``): the sum of the areas of their stem elements (stem organs) or leaf
elements (blades).


Leaf shapes (``leaf_shapes``)
=============================

**The width profile** of a class is linear between its points :math:`(s_i, r_i)`, sorted by :math:`s`. The area of a
part (``LeafShapes.relative_area``) is the exact integral of the linear profile, from the cumulative integrals at the
points computed once (``LeafShapes._cumulative``): for :math:`x` in :math:`[s_i, s_{i+1}]`,

.. math::

   \int_0^x r = C_i + r_i\,(x - s_i) + \tfrac{1}{2}\,\frac{r_{i+1} - r_i}{s_{i+1} - s_i}\,(x - s_i)^2, \qquad
   C_i = \sum_{j < i} \tfrac{1}{2}\,(r_j + r_{j+1})(s_{j+1} - s_j),

the bounds clipped to :math:`[0, 1]` and to the profile's points. **The mean width** of a part
:math:`[s_b, s_t]` (``LeafShapes.mean_width``, Adel's ``blade_elt_width``) is the mean of :math:`r` at its two bounds
and at the profile's :math:`m` points strictly inside it,

.. math::

   \bar r = \frac{r(s_b) + r(s_t) + \sum_{s_b < s_i < s_t} r_i}{2 + m},

a point within :math:`10^{-9}` of a bound counted as the bound (``BOUND_TOLERANCE``).

**Shapes.** A class is numbered from 1 in the database's classes' ascending numerical order (Adel's ``Lindex``); a
shape in its class from 1 (``Lseed``), among the shapes that could be fitted; a dynamic database (echap's) has one
fitted midrib per age class. ``LeafShapes.shape_index`` gives each blade its fitted midrib (−1: none). **The age
class** (``LeafShapes.age_class``, Adel's ``get_age_index``) of a blade of age :math:`a` on an axis of phyllochron
:math:`\varphi` (``phyllochron`` for the main stem, ``tiller_phyllochron`` for every tiller) is the bin of
:math:`a/\varphi - 0.3` among the bins :math:`\beta_0 < \dots < \beta_m`: the index :math:`j` with
:math:`\beta_j < a/\varphi - 0.3 \le \beta_{j+1}`, 0 below the bins, :math:`m - 1` above them; always 0 for a static
database. echap's bins (``LeafShapes.echap``) lie half-way between the measured ages, the first one below the first
age minus 1, the last one above the last age plus 1.

**Fitting** (``leaf_shapes.fit_midrib``, at the database's construction), as Adel's ``fitting.fit_leaves``:

1. the midrib's points :math:`(x, y)` by a smoothing spline (``splprep``, smoothing 0.01, order 3; failing that
   scipy's defaults, then order 1), evaluated at 100 parameters, its curvilinear abscissa :math:`s`; :math:`x`,
   :math:`y`, :math:`s` divided by the midrib's length (``_smooth``);
2. the width profile by a smoothing spline (smoothing 0.0001, order 3; failing that scipy's defaults) at 200
   parameters, :math:`s` divided by its maximum, :math:`r` set to :math:`\max(r, 0)` divided by its maximum, then
   interpolated at the midrib's abscissae;
3. the polyline of points :math:`(x, r, y)` simplified (``_simplification_kept``, Adel's ``simplification.cost``):
   the interior point nearest the line through its two neighbours removed first, by the squared distance
   :math:`|(\mathbf{p} - \mathbf{a}) \times (\mathbf{b} - \mathbf{a})|^2 / |\mathbf{b} - \mathbf{a}|^2`, its
   neighbours' distances updated, until ``discretisation_level`` interior points remain (9 for the samples, 7 for
   echap's); then the second and the last but one points removed;
4. the kept points' abscissa, :math:`x` and :math:`y` multiplied by :math:`\max s_\text{smooth} / \max s_\text{kept}`,
   and :math:`r` by the ratio of the smooth shape's area (Simpson's rule) to the polygon's (trapezoids)
   (``_simplify``);
5. the tip set to :math:`s = 1`, :math:`r = 0`. A shape that cannot be fitted is left out, the next shapes of its
   class taking its number; a dynamic shape is left out when one of its age classes cannot be fitted.


The placement (``turtle``)
==========================

**The frame operations** (metafspm's ``coupling.geometry``, PlantGL's turtle): a frame is a heading :math:`H` and an
up :math:`U`, its left :math:`L = U \times H`; :math:`R_{\mathbf{a}}(\theta)` is the rotation by :math:`\theta` about
the unit axis :math:`\mathbf{a}` (Rodrigues, right-handed).

.. math::

   \text{rollR}(\theta):\ U \leftarrow R_H(-\theta)\,U; \qquad
   \text{down}(\theta):\ H \leftarrow R_L(\theta)\,H,\ U \leftarrow R_L(\theta)\,U; \qquad
   \text{rollToVert}:\ U \leftarrow H \times \frac{(-H_y, H_x, 0)}{\|(-H_y, H_x, 0)\|},

with :math:`(0, -1, 0)` for the normalised vector when :math:`H` is vertical (PlantGL's convention).

**The traversal** (``turtle.traversal_order``, once per topology): plant by plant from its base element, depth first,
a node's branches before its successors. A tiller borne on the main stem's phytomer 1 (T1) is visited after that
phytomer's sheath and blade, as the last branch of its internode's last element; its first element takes the heading
and up of the element it is branched on (``frame_source``). It then starts at the coleoptile's top, where phytomer 2
and T2 start. The MPG keeps T1 on phytomer 1.

**The turtle** (``turtle.interpret``, compiled; Adel's ``AdelVisitor``). It does not follow the topology: each element
starts from the frame on top of the elements placed so far (the top), with two saved frames, the main stem's top and
the tillers' base:

- a plant's base element: the top at the plant's position :math:`(x, y, z)`, :math:`H = (0, 0, 1)`,
  :math:`U = (-1, 0, 0)`, then rollR(rotation); both saved frames set to it;
- the first element of a tiller's first phytomer, at the axis' first visit: after a main stem element, the main
  stem's top and the tillers' base are saved as the top (with T1's heading and up from its ``frame_source``); after
  another axis' element, the top returns to the tillers' base;
- any other change of axis: the top returns to the main stem's top saved (a tiller's continuation after its own
  tiller: the frame saved at the tiller's start, with T1's heading and up from its ``frame_source``, as the legacy
  saves the main stem's top before the coleoptile);
- with ``own_axis_top`` 1 (correction AD-37), at a change of axis: an axis visited before returns to its own top
  where the turtle left it (a tiller with the flag of its first inclined stem element as it left it), a tiller's
  first element starts on the top its topological parent (its bearer's element) left;
- a stem element with an inclination :math:`\theta \ne 0`: the axis' first one, on a tiller, is rolled by the axis'
  azimuth (rollR(``axis_azimuth``)), saved as the tillers' base, then turned down(:math:`\theta`); the next ones are
  rolled to the vertical (rollToVert), turned down(:math:`\theta`), then rolled back into their azimuth plane by the
  angle :math:`\gamma` between their up before and after rollToVert: rollR(:math:`\gamma`) when rollToVert raised
  the left's :math:`z`, rollR(:math:`-\gamma`) otherwise;
- a stem element of an internode with an azimuth :math:`\alpha \ne 0` (``organ_azimuth``): rollR(:math:`\alpha`);
- each element's frame (``x1``, ``y1``, ``z1``, ``heading_*``, ``up_*``) is the frame so obtained; then a stem
  element moves the top forward by its length, a blade element by its rolled part when it has one; hidden elements
  leave the top as it is. Its end (``x2``, ``y2``, ``z2``) is its position plus its length (a blade element: its
  rolled part) along :math:`H`.

The inclinations are 0 unless ``tiller_inclination`` is 1 (AdelDyn's vertical axes). With it, static Adel's rules
(Adel.R) give the axes' first inclinations and azimuths (``ShootArchitecture._first_inclination``, an axis' step) and
the stem organs' inclinations (``ShootArchitecture._inclination``, an organ's step), each step getting its parameters
per axis or per organ, its plant's:

**The base inclination** (``turtle.base_inclination``, Adel.R's ``getIncT``): an axis' first internode is inclined by

.. math::

   \theta_0 = \begin{cases} \text{incT} & \text{main stem} \\
   3° & \text{a tiller, } HS \le HS_0 \\
   \max\big(3°,\ \min(\text{incT},\ \kappa\,(HS - HS_0))\big) & \text{a tiller, } HS > HS_0 \end{cases}

with the axis' Haun stage :math:`HS` (``haun_stage``, an input), :math:`HS_0` ``inclination_start_haun_stage`` (1)
and :math:`\kappa` ``inclination_rate`` (30° per Haun stage): the axis' ``first_inclination``.

**The straightening** (``turtle.stem_inclinations``, Adel.R's ``axe_inclination``): per axis, its stem organs sorted
by rank, the internode before the sheath, those of positive visible length :math:`\ell_j` in turn, with
:math:`\alpha = 90° - \theta_0`, the cumulated length :math:`\Lambda_j = \sum_{i \le j} \ell_i` and the horizontal
distance :math:`D_j = \Lambda_j \cos\alpha`. At the first :math:`j` with :math:`D_j \ge D` (:math:`D` the
straightening distance, dredT; none when :math:`D = 0`):

.. math::

   j > 0:\quad \beta = \arccos\frac{D - D_{j-1}}{\Lambda_j - \Lambda_{j-1}}, \quad \theta_j = -(\beta - \alpha), \quad
   \theta_{j+1} = -(90° - \beta);

.. math::

   j = 0:\quad \beta = \arccos\frac{D}{\Lambda_0}, \quad \theta_0 = 90° - \beta, \quad \theta_1 = -\theta_0 .

The other organs' inclinations are 0, the first internode's :math:`\theta_0`. **The azimuths:** each axis' first
internode's ``organ_azimuth`` becomes its ``first_azimuth`` (azT), each axis' ``axis_azimuth`` its
``insertion_azimuth`` (azTb).


The triangles (``ShootArchitecture.triangles``, ``blade_mesh``)
===============================================================

Each element is drawn in its frame: a local point :math:`(p_x, p_y, p_z)` is placed at
:math:`\mathbf{x} + p_x U + p_y\,(H \times U) + p_z H` (metafspm's ``coupling.geometry.place``).

**Stem elements** of length :math:`\ell \ge` ``min_length`` (Adel's 0.01 cm, 1e-4 m) and **blades' rolled parts**
(:math:`\rho_k > 0`): a slim cylinder (``coupling.geometry.slim_cylinders``, Adel's ``slim_cylinder``) of radius
:math:`d/2` (the rolled diameter's half for a rolled part): three points at the angles :math:`0, 2\pi/3, 4\pi/3` at
its base, three at the same angles plus :math:`\pi` at :math:`p_z = \ell`, eight triangles.

**Blade elements' flat parts** (``ShootArchitecture._blade_triangles``, ``blade_mesh.blade_element_triangles``), for
a blade of visible length at least ``min_length`` and an element with a flat part, from its shape's fitted midrib
:math:`(x_i, y_i, s_i, r_i)` (unit length):

1. **inclination** (``blade_mesh.incline``, Adel's ``incline_leaf``): with the base angle
   :math:`\phi = \operatorname{atan2}(x_1 - x_0,\ y_1 - y_0)` and a factor :math:`c`, the midrib rotated by
   :math:`\psi = \phi - \min(\pi, c\,\phi)`,

   .. math::

      x' = x_0 + x\cos\psi - y\sin\psi, \qquad y' = y_0 + x\sin\psi + y\cos\psi,

   :math:`c = 1` for a dynamic database (its shapes carry their inclination), the blade's ``inclination`` otherwise;
   then :math:`x' \leftarrow -x'`;
2. **the part** (``blade_mesh.blade_part``, Adel's ``fitting.leaf_element``): with the visible fraction
   :math:`v = \min(l_v, L)/L` and the part :math:`[\sigma_b, \sigma_t]` = [flat_base, flat_top]: the whole shape
   scaled when :math:`v = 1` and the part is :math:`[0, 1]`; otherwise :math:`r` at the tip set to 0, the midrib's
   abscissae the shape's points below :math:`v` and those above :math:`1 - v` moved down by :math:`1 - v`, with
   :math:`v\sigma_b` and :math:`v\sigma_t`, kept within :math:`[v\sigma_b, v\sigma_t]`; a point nearer than
   :math:`v(\sigma_t - \sigma_b)/(2 n_s)` to the previous one dropped (:math:`n_s` the shape's number of points), the
   two ends kept; at each abscissa :math:`u`, the midrib point :math:`(x'(u), y'(u))` and the width
   :math:`r(u + 1 - v)` (the blade grows out of the tube: its visible top is the shape's top), interpolated linearly;
   the points cut after the first zero width; :math:`x, y` multiplied by :math:`L`, :math:`r` by :math:`W`;
3. **the strip** (``blade_mesh.strip_triangles``, Adel's ``leaf_to_mesh`` without twist): the points
   :math:`(x_i, -r_i/2, y_i)` and :math:`(x_i, r_i/2, y_i)`, :math:`i = 0 \dots m - 1`, two triangles between
   consecutive points up to :math:`m - 2`, one at the tip when its width is below 1e-5 m (Adel's 0.001 cm), two
   otherwise;
4. the triangles of area at most ``min_length``\ :sup:`2` / 100 left out; the mesh lifted by the element's rolled
   part along :math:`p_z`, above its rolled cylinder.

The blade elements of one blade are drawn from frames that only their rolled parts move, each its own part of the
whole visible midrib. An element's local mesh is kept between calls, by MPG vertex, while what it is made from is the
same (its shape, inclination factor, mature length and width, visible length, flat part, rolled part, the area
threshold); every element is placed again at each call.


The scenario development (``development.ScenarioDevelopment``)
==============================================================

A deterministic scenario, linear in thermal time, not Adel-Wheat's kinetics. A plant's thermal time :math:`t` starts
at ``start`` (by default its main stem's ``emf1``) minus its ``delay``, and is the plant's at the start of each call
(``ScenarioDevelopment._degree_days``): with the scene's forcing ``air_temperature`` :math:`T` (°C), a call adds the
degree-days of the step before it, :math:`\max(0, T - T_b)\,\Delta t / 1\,\text{d}`, :math:`T` read at that step's
start, :math:`T_b` ``base_temperature`` (0 °C), :math:`\Delta t` the time step; without it, a call adds
``thermal_time_step``. With the phyllochron :math:`\varphi` (``phyllochron``, every axis), an axis' first leaf at :math:`e`
(``emf1``), rank :math:`n` appearing at :math:`a_n = e + (n - 1)\varphi`, and
:math:`\operatorname{ramp}(b, D) = \operatorname{clip}((t - b)/D, 0, 1)`:

.. math::

   \texttt{axis\_number} = \#\{\text{axes}: t \ge e\}, \qquad
   \texttt{phytomer\_number} = \#\{n \le n_f: t \ge a_n\}, \qquad
   HS = \operatorname{clip}\big((t - e)/\varphi,\ 0,\ n_f\big),

.. math::

   l_\text{blade} &= \text{Ll}_n\,\operatorname{ramp}(a_n,\ 1.6\varphi), \\
   l_\text{sheath} &= \text{Gl}_n\,\operatorname{ramp}(a_n + 1\varphi,\ 1\varphi), \\
   l_\text{internode} &= \text{El}_n\,\operatorname{ramp}(a_n + 2\varphi,\ 1.5\varphi),

the durations' factors being the parameters ``blade_duration``, ``sheath_start``, ``sheath_duration``,
``internode_start``, ``internode_duration`` (their defaults shown). The previous sheath is
:math:`l'_\text{sheath} = \text{Gl}_{n-1}\,\operatorname{ramp}(a_{n-1} + 1\varphi, 1\varphi)`; the visible lengths
are :math:`\max(0, l_\text{blade} - l'_\text{sheath})`, :math:`\max(0, l_\text{sheath} - l'_\text{sheath})` and
:math:`0.5\,l_\text{internode}` (``internode_visible_fraction``); the blade's senesced length

.. math::

   l_s = \min\big(\text{Ll}_n\,\operatorname{ramp}(a_n + 5\varphi,\ 3\varphi)\,[l_\text{blade} > 0],\ l_\text{blade}\big)

(``senescence_start``, ``senescence_duration``), and its age :math:`t - a_n` once it has appeared. The lengths are
multiplied by the tables' unit. The scenario writes no rolled length. Each output location has its steps, which get
their parameters there, a plant's value per axis or organ: per plant the thermal time and ``axis_number``
(``_thermal_time``, ``_axis_number``), per axis ``phytomer_number`` and :math:`HS` (``_axis_development``), per organ
the lengths and the age (``_organ_development``). ``ScenarioShoot`` runs the development and the shoot together, in
one schedule (metafspm's ``CompositeModel.run_together``): the development's potential steps, the shoot's
segmentation (the axes and phytomers asked for), the development's postsegmentation steps (the new axes' Haun stages
and the new organs' lengths, ``ScenarioDevelopment._new_axes``, ``_new_organs``), then the shoot's
postsegmentation.
