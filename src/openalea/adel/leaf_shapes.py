"""
Leaf shape databases for ShootArchitecture: each shape class' width profile r(s) along the blade (s from 0 at the base
to 1 at the tip, r relative to the maximal width), and the areas of blade parts, vectorised.

    shapes = LeafShapes.samples()                        Adel's samples (data/SRSo.csv)
    shapes = LeafShapes.echap(sr_model="Mercia")         echap's (echap_leaf_data), the ones cnwgrass uses
    shapes.relative_area(classes, s_base, s_top)         the integral of r between s_base and s_top, per blade part

A shape class is numbered from 1, in the database's classes' ascending order (Adel's Lindex). The width profile is
linear between its points (as every Adel function interpolates it), and areas are its exact integral: cumulative at
the points once, so that any part's area is a difference of two interpolated values.
"""
import os

import numpy as np
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))


class LeafShapes:
    """Width profiles by shape class, {class id: (s, r)}, with their cumulative integrals."""

    def __init__(self, profiles):
        self.classes = sorted(profiles)                     # class index k (from 1) is self.classes[k - 1]
        self.profiles = {key: (np.asarray(s, dtype=float), np.asarray(r, dtype=float))
                         for key, (s, r) in profiles.items()}
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
    def from_table(cls, table):
        """From a table of columns rankclass, s, r."""
        return cls({int(key): (group["s"].to_numpy(), group["r"].to_numpy())
                    for key, group in table.groupby("rankclass")})

    @classmethod
    def samples(cls):
        """Adel's sample profiles (data/SRSo.csv: 6 classes)."""
        return cls.from_table(pd.read_csv(os.path.join(DATA, "data", "SRSo.csv")))

    @classmethod
    def echap(cls, sr_model="Mercia"):
        """echap's profiles (echap_leaf_data/srdb_Grignon*.csv: 2 classes, the base leaves and the top ones)."""
        files = {"Mercia": "srdb_GrignonMercia2010.csv", "Tremie": "srdb_GrignonTremie2011.csv"}
        table = pd.read_csv(os.path.join(DATA, "echap_leaf_data", files[sr_model]), names=["rankclass", "s", "r"],
                            skiprows=1)
        return cls.from_table(table)

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
        first, stop = self._search(classes, s_base, "right"), self._search(classes, s_top, "left")
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
