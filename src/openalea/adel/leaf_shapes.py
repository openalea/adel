"""
Leaf shape databases for ShootArchitecture: each shape class' width profile r(s) along the blade (s from 0 at the base
to 1 at the tip, r relative to the maximal width), the areas of blade parts, vectorised; and the blades' midribs,
fitted once for their meshes.

    shapes = LeafShapes.samples()                        Adel's samples (data/SRSo.csv, data/So99.csv)
    shapes = LeafShapes.echap(sr_model="Mercia")         echap's (echap_leaf_data), the ones cnwgrass uses
    shapes.relative_area(classes, s_base, s_top)         the integral of r between s_base and s_top, per blade part
    shapes.shape_index(classes, seeds, age_classes)      the fitted midrib of each blade (-1: none)
    shapes.fitted[i]                                     a fitted midrib: (x, y, s, r) arrays

A shape class is numbered from 1, in the database's classes' ascending order (Adel's Lindex); a shape in its class
from 1 (Adel's Lseed); a dynamic database (echap's) has a midrib per age class, the class of a blade's age in
phyllochrons (age_class). The width profile is linear between its points (as every Adel function interpolates it),
and areas are its exact integral: cumulative at the points once, so that any part's area is a difference of two
interpolated values.

The midribs are fitted as Adel fits them (fitting.fit_leaves: fit2, then simplify to discretisation_level points): a
smoothing spline of the midrib, resampled and scaled to a unit length, the width profile's spline on its curvilinear
abscissa, then the polyline simplified (the points removed one by one, the one nearest the line between its
neighbours first), its width rescaled to keep the smooth shape's area.
"""
import heapq
import os

import numpy as np
import pandas as pd
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import splev, splprep

DATA = os.path.dirname(os.path.abspath(__file__))
BOUND_TOLERANCE = 1e-9                  # relative abscissa


class LeafShapes:
    """Width profiles by shape class, {class id: (s, r)}, with their cumulative integrals; midribs by class, {class
    id: [per shape: (x, y)]}, or for a dynamic database [per shape: {age class: (x, y)}] with the age classes' bounds
    *age_bins* (in phyllochrons)."""

    def __init__(self, profiles, midribs=None, age_bins=None, discretisation_level=9):
        self.classes = sorted(profiles)                     # class index k (from 1) is self.classes[k - 1]
        self.profiles = {key: (np.asarray(s, dtype=float), np.asarray(r, dtype=float))
                         for key, (s, r) in profiles.items()}
        self.dynamic = age_bins is not None
        self.age_bins = None if age_bins is None else np.asarray(age_bins, dtype=float)
        self.discretisation_level = discretisation_level
        self._fit_midribs(midribs or {})
        # All the profiles' points one after the other: per class, its points' positions and cumulative integrals
        s_all, r_all, c_all, offsets = [], [], [], [0]
        for key in self.classes:
            s, r = self.profiles[key]
            order = np.argsort(s, kind="stable")
            s, r = s[order], r[order]
            c = np.concatenate([[0.], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(s))])
            s_all.append(s)
            r_all.append(r)
            c_all.append(c)
            offsets.append(offsets[-1] + s.size)
        self._s, self._r, self._c = (np.concatenate(a) for a in (s_all, r_all, c_all))
        self._offsets = np.asarray(offsets)
        self._r_sums = np.concatenate([[0.], np.cumsum(self._r)])     # the points' r summed, for mean widths

    @classmethod
    def from_table(cls, table, midribs=None, age_bins=None, discretisation_level=9):
        """From a table of columns rankclass, s, r (and the midribs, as LeafShapes takes them)."""
        return cls({int(key): (group["s"].to_numpy(), group["r"].to_numpy())
                    for key, group in table.groupby("rankclass")}, midribs, age_bins, discretisation_level)

    @classmethod
    def samples(cls, discretisation_level=9):
        """Adel's sample profiles (data/SRSo.csv: 6 classes) and midribs (data/So99.csv: shapes by class, its
        column rank, in their column lindex' order)."""
        table = pd.read_csv(os.path.join(DATA, "data", "So99.csv"))
        midribs = {}
        for (rank, _), group in table.groupby(["rank", "lindex"]):
            midribs.setdefault(int(rank), []).append((group["x"].to_numpy(), group["y"].to_numpy()))
        return cls.from_table(pd.read_csv(os.path.join(DATA, "data", "SRSo.csv")), midribs,
                              discretisation_level=discretisation_level)

    @classmethod
    def echap(cls, sr_model="Mercia", xy_model="Tremie_byleafclass", discretisation_level=7):
        """echap's profiles (echap_leaf_data/srdb_Grignon*.csv: 2 classes, the base leaves and the top ones) and its
        median midribs by age (echap_leaf_data/median_leaf_trajectories_*.csv: one shape per class, an age class
        per age measured, bounded half-way between them)."""
        files = {"Mercia": "srdb_GrignonMercia2010.csv", "Tremie": "srdb_GrignonTremie2011.csv"}
        table = pd.read_csv(os.path.join(DATA, "echap_leaf_data", files[sr_model]), names=["rankclass", "s", "r"],
                            skiprows=1)
        trajectories = {"MerciaRht_byleafclass": "MerciaRht_Grignon2010_byleafclass",
                        "Tremie_byleafclass": "Tremie_Boigneville2012_2013_byleafclass",
                        "Soissons_byleafclass": "Soissons_Grignon2010_byleafclass",
                        "MerciaRht": "MerciaRht_Grignon2010",
                        "Tremie": "Tremie_Boigneville2012_2013", "Soissons": "Soissons_Grignon2010",
                        "Tremie12": "Tremie12", "Tremie13": "Tremie13"}
        points = pd.read_csv(os.path.join(DATA, "echap_leaf_data",
                                          f"median_leaf_trajectories_{trajectories[xy_model]}.csv"),
                             names=["age", "x", "y", "lindex"], skiprows=1)
        ages = np.array(sorted(set(points["age"])))
        middles = (ages[1:] + ages[:-1]) / 2.
        age_bins = np.concatenate([[ages[0] - 1.], middles, [ages[-1] + 1.]])
        points["age_class"] = pd.cut(points["age"], age_bins, labels=False)
        midribs = {}
        for key in sorted(set(points["lindex"])):
            mine = points.loc[points["lindex"] == key]
            midribs[int(key)] = [{int(age): (group["x"].to_numpy(), group["y"].to_numpy())
                                  for age, group in mine.groupby("age_class")}]
        return cls.from_table(table, midribs, age_bins, discretisation_level)

    # --- The midribs ---

    def _fit_midribs(self, midribs):
        """Each shape's fitted midrib, in self.fitted, and the index of each (class, shape, age class) in it (the
        shapes that cannot be fitted left out, the next ones of their class taking their numbers, as Adel's
        fit_leaves does; a dynamic shape is left out when one of its age classes cannot be fitted)."""
        self.fitted = []
        n_ages = 1 if not self.dynamic else self.age_bins.size - 1
        kept = {}
        for key in self.classes:
            s, r = (np.asarray(values, dtype=float) for values in self.profiles[key])     # in the file's order
            shapes = []
            for shape in midribs.get(key, []):
                by_age = shape if self.dynamic else {0: shape}
                fitted = {age: fit_midrib(x, y, s, r, self.discretisation_level) for age, (x, y) in by_age.items()}
                if all(f is not None for f in fitted.values()):
                    shapes.append(fitted)
            kept[key] = shapes
        n_seeds = max([len(shapes) for shapes in kept.values()] + [0])
        self._shape_table = np.full((len(self.classes) + 1, n_seeds + 1, n_ages), -1, dtype=np.int64)
        for k, key in enumerate(self.classes, start=1):
            for seed, fitted in enumerate(kept[key], start=1):
                for age, shape in fitted.items():
                    self._shape_table[k, seed, age] = len(self.fitted)
                    self.fitted.append(shape)

    def age_class(self, age_in_phyllochrons):
        """A blade's age class from its age in phyllochrons (Adel's get_age_index): the bin below it, the first one
        below the bins, the last one above them; 0 for a static database."""
        age = np.asarray(age_in_phyllochrons, dtype=float)
        if not self.dynamic:
            return np.zeros(age.shape, dtype=np.int64)
        below = np.maximum(np.searchsorted(self.age_bins, age) - 1, 0)
        return np.where(below >= self.age_bins.size - 1, below - 1, below)

    def shape_index(self, classes, seeds, age_classes):
        """Per blade, its fitted midrib's index in self.fitted (-1: none, or no such shape)."""
        classes, seeds = np.asarray(classes, dtype=np.int64), np.asarray(seeds, dtype=np.int64)
        age_classes = np.asarray(age_classes, dtype=np.int64)
        table = self._shape_table
        valid = ((classes >= 1) & (classes < table.shape[0]) & (seeds >= 1) & (seeds < table.shape[1])
                 & (age_classes >= 0) & (age_classes < table.shape[2]))
        out = np.full(classes.shape, -1, dtype=np.int64)
        out[valid] = table[classes[valid], seeds[valid], age_classes[valid]]
        return out

    def _search(self, classes, x, side):
        """Per point, the index (in the concatenated points) where x would be inserted in its class' points."""
        out = np.empty(np.shape(x), dtype=np.int64)
        for k in np.unique(classes):
            these = classes == k
            s = self._s[self._offsets[k - 1]:self._offsets[k]]
            out[these] = self._offsets[k - 1] + np.searchsorted(s, x[these], side=side)
        return out

    def _interpolate(self, classes, x):
        """r at x, per point."""
        out = np.empty(np.shape(x))
        for k in np.unique(classes):
            these = classes == k
            span = slice(self._offsets[k - 1], self._offsets[k])
            out[these] = np.interp(x[these], self._s[span], self._r[span])
        return out

    def mean_width(self, classes, s_base, s_top):
        """Adel's width of a blade part, relative to the maximal width (geometric_elements.blade_elt_width): the
        mean of r at s_base, s_top and the profile's points strictly between them."""
        classes = np.asarray(classes, dtype=np.int64)
        s_base = np.clip(np.asarray(s_base, dtype=float), 0., 1.)
        s_top = np.clip(np.maximum(np.asarray(s_top, dtype=float), s_base), 0., 1.)
        if s_base.size == 0:
            return np.zeros(0)
        # strictly inside: a point within BOUND_TOLERANCE of a bound is the bound itself (the bounds of sectors fall
        # on the profile's points, and which side rounding puts them on would decide)
        first = self._search(classes, s_base + BOUND_TOLERANCE, "right")
        stop = self._search(classes, s_top - BOUND_TOLERANCE, "left")
        inside = np.maximum(stop - first, 0)
        total = self._r_sums[np.maximum(stop, first)] - self._r_sums[first]
        return (self._interpolate(classes, s_base) + self._interpolate(classes, s_top) + total) / (2 + inside)

    def _cumulative(self, classes, x):
        """The integral of r from 0 to x, per point (classes from 1)."""
        classes = np.asarray(classes, dtype=np.int64)
        x = np.asarray(x, dtype=float)
        start, stop = self._offsets[classes - 1], self._offsets[classes]
        s_first, s_last = self._s[start], self._s[stop - 1]
        x = np.clip(x, s_first, s_last)
        # the segment of each point: its class' points, searched in the class' slice
        position = np.empty(x.shape, dtype=np.int64)
        for k in np.unique(classes):
            these = classes == k
            s = self._s[self._offsets[k - 1]:self._offsets[k]]
            position[these] = self._offsets[k - 1] + np.clip(np.searchsorted(s, x[these], side="right") - 1, 0,
                                                             s.size - 2)
        s0, s1, r0, r1 = self._s[position], self._s[position + 1], self._r[position], self._r[position + 1]
        slope = np.where(s1 > s0, (r1 - r0) / np.where(s1 > s0, s1 - s0, 1.), 0.)
        dx = x - s0
        return self._c[position] + r0 * dx + 0.5 * slope * dx ** 2

    def relative_area(self, classes, s_base, s_top):
        """The integral of r between s_base and s_top (clipped to [0, 1], s_top at least s_base), per blade part;
        times the blade's mature length and maximal width, its area."""
        s_base = np.clip(np.asarray(s_base, dtype=float), 0., 1.)
        s_top = np.clip(np.maximum(np.asarray(s_top, dtype=float), s_base), 0., 1.)
        if np.size(s_base) == 0:
            return np.zeros(0)
        return self._cumulative(classes, s_top) - self._cumulative(classes, s_base)


def curvilinear_abscissa(x, y):
    """The cumulative length along the polyline (x, y), from 0."""
    out = np.zeros(len(x))
    out[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    return out.cumsum()


def fit_midrib(x, y, s, r, number_of_points):
    """A midrib fitted as Adel's fit_leaves fits it (fitting.fit3, then its tip at s, r = 1, 0): (x, y, s, r)
    arrays, the midrib of length 1 (its curvilinear abscissa s), r relative to the maximal width; None when it
    cannot be fitted."""
    try:
        shape = _simplify(*_smooth(np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                   np.asarray(s, dtype=float), np.asarray(r, dtype=float)), number_of_points)
    except Exception:                                   # as Adel's: a shape that cannot be fitted is left out
        return None
    shape[2][-1], shape[3][-1] = 1., 0.
    return shape


def _smooth(x, y, s, r):
    """Adel's fitting.fit2: the midrib's smoothing spline at 100 points, scaled to a unit length; the width profile's
    spline at 200 points (r positive, relative to its maximum), interpolated on the midrib's abscissa."""
    smooth, order = 0.01, 3
    attempts = (dict(s=smooth, k=order, nest=-1), dict(), dict(k=1))     # Adel's, in its order
    for attempt, options in enumerate(attempts):
        try:
            midrib, _ = splprep([x, y], **options)
            break
        except Exception:
            if attempt == len(attempts) - 1:
                raise
    x_new, y_new = splev(np.linspace(0, 1, 100), midrib)
    s_new = curvilinear_abscissa(x_new, y_new)
    length = s_new.max()
    s_new /= length
    x_new /= length
    y_new /= length
    try:
        profile, _ = splprep([s, r], s=smooth / 100., k=order, nest=-1)
    except Exception:
        profile, _ = splprep([s, r])
    s_profile, r_profile = splev(np.linspace(0, 1, 200), profile)
    s_profile = s_profile / s_profile.max()
    r_profile = np.maximum(r_profile, 0) / r_profile.max()
    return x_new, y_new, s_new, np.interp(s_new, s_profile, r_profile)


def _squared_distance(point, start, end):
    """The squared distance of *point* to the line through *start* and *end*."""
    line = end - start
    return float(np.sum(np.cross(point - start, line) ** 2) / np.sum(line ** 2))


def _simplification_kept(points, number_of_points):
    """The points Adel's simplification.cost keeps: the interior points removed one at a time, the one nearest the line
    through its neighbours first, until number_of_points remain inside; then, as Adel's, the second and the last but
    one removed too (a mask)."""
    n = len(points)
    kept = np.ones(n, dtype=bool)
    neighbours = [[i - 1, i + 1] for i in range(n)]
    heap, cost = [], {}
    for i in range(1, n - 1):
        cost[i] = _squared_distance(points[i], points[i - 1], points[i + 1])
        heapq.heappush(heap, [cost[i], i])
    while len(heap) > number_of_points:
        _, i = heapq.heappop(heap)
        kept[i] = False
        left, right = neighbours[i]
        for side, other in ((left, right), (right, left)):
            if side in (0, n - 1):
                continue
            neighbours[side][1 if side == left else 0] = other
            beyond = neighbours[side][0 if side == left else 1]
            start, end = (beyond, other) if side == left else (other, beyond)
            distance = _squared_distance(points[side], points[start], points[end])
            del heap[heap.index([cost[side], side])]
            cost[side] = distance
            heapq.heappush(heap, [distance, side])
    kept[1] = kept[-2] = False
    return kept


def _simplify(x, y, s, r, number_of_points):
    """Adel's fitting.simplify: the polyline (x, r, y) simplified to number_of_points points inside, rescaled to the
    smooth shape's abscissa, its width rescaled to the smooth shape's area (Simpson's rule on the smooth one, the
    trapezoids on the polygon)."""
    kept = _simplification_kept(np.stack([x, r, y], axis=1), number_of_points)
    x_kept, y_kept, r_kept = x[kept], y[kept], r[kept]
    s_kept = curvilinear_abscissa(x_kept, y_kept)
    adjust = max(s) / max(s_kept)
    s_kept, x_kept, y_kept = s_kept * adjust, x_kept * adjust, y_kept * adjust
    r_kept = r_kept * (simpson(r, x=s) / trapezoid(r_kept, x=s_kept))
    return [x_kept, y_kept, s_kept, r_kept]
