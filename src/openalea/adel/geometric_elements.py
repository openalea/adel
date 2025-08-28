"""
New proposal for computing organ shapes
"""

import numpy
import pandas
import os
from scipy.integrate import simpson

import openalea.plantgl.all as pgl
from math import radians, pi, cos, sin
import openalea.adel.fitting as fitting
from functools import reduce
from copy import deepcopy
datadir = os.path.dirname(__file__)


def genGeoLeaf(nlim=4, dazt=60, dazb=10):
    """generate geoLeaf function for Adel"""
    rcode = """
    geoLeaf <- list(
     Azim = function(a,n,nf) {{
            ntop = nf - n
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,nf) {{
              nf - n + 1}}
              )
        """
    return rcode.format(ntoplim=nlim, dazTop=dazt, dazBase=dazb)


def leaf_keys(lindex, lseed, db):
    """convert R-style lindex/lseed (also called LcType/Lindex in canopy table)
    into (keys,index) of python xy/sr data bases
    """
    if 1 > lindex or lindex > len(db) or lseed < 1:
        raise KeyError("invalid index for leaf shape database")
    keys = sorted(db.keys())
    return keys[lindex - 1], lseed - 1


def xydb_to_csv(xydb, filename):
    dat = [
        (numpy.repeat(k, len(x)), numpy.repeat(i, len(x)), x, y)
        for k in xydb
        for i, (x, y) in enumerate(xydb[k])
    ]
    dat = list(reduce(lambda x, y: x + y, map(lambda x: zip(*x), dat), []))
    df = pandas.DataFrame.from_records(dat, columns=("rank", "lindex", "x", "y"))
    df.to_csv(filename, index=False, sep=",", decimal=".")


def xydb_from_csv(filename):
    df = pandas.read_csv(filename, sep=",", decimal=".")
    grouped = df.groupby(["rank", "lindex"])
    xydb = {str(int(r)): [] for r in set(df["rank"])}
    for (rank, lindex), data in iter(grouped):
        xydb[str(int(rank))] += [
            (numpy.array(data.loc[:, "x"]), numpy.array(data.loc[:, "y"]))
        ]
    return xydb


def srdb_to_csv(srdb, filename):
    dat = [(numpy.repeat(k, len(x)), x, y) for k, (x, y) in srdb.iteritems()]
    dat = reduce(lambda x, y: x + y, map(lambda x: zip(*x), dat), [])
    df = pandas.DataFrame.from_records(dat, columns=("rankclass", "s", "r"))
    df.to_csv(filename, index=False, sep=",", decimal=".")


def srdb_from_csv(filename):
    df = pandas.read_csv(filename, sep=",", decimal=".")
    grouped = df.groupby("rankclass")
    srdb = {
        str(int(r)): (numpy.array(data.loc[:, "s"]), numpy.array(data.loc[:, "r"]))
        for r, data in grouped
    }
    return srdb


def curvilinear_abscisse(x, y):
    s = numpy.zeros(len(x))
    s[1:] = numpy.sqrt(numpy.diff(x) ** 2 + numpy.diff(y) ** 2)
    return s.cumsum()


# function for handling dynamical leaf shape databases
def incline_leaf(shape, inclin, relative_angle=True):
    """transform a xysr tuple representing leaf shape to get a given angle at leaf base.
    - angle the desired angle (deg)
    - if relative_angle == True, angle is interpreted as a multiplier to original shape angle
    """
    Linc = inclin
    x, y = shape[0], shape[1]
    init_angle = pgl.angle((x[1] - x[0], y[1] - y[0]), (0, 1))

    if relative_angle:
        angle = Linc * init_angle
        angle = min(pi, angle)
    else:
        angle = radians(Linc)

    rotation_angle = init_angle - angle

    # rotation of the midrib
    cos_a = cos(rotation_angle)
    sin_a = sin(rotation_angle)

    x1 = x[0] + cos_a * x - sin_a * y
    y1 = y[0] + sin_a * x + cos_a * y
    leaf = x1, y1, shape[2], shape[3]

    return leaf


class Leaves:
    def __init__(
        self,
        xydb=None,
        srdb=None,
        geoLeaf=None,
        dynamic_bins=None,
        discretisation_level=9,
        twist=0,
    ):
        self.leaves = None
        if xydb is None:
            data = datadir + "/data/So99.csv"
            xydb = xydb_from_csv(data)

        if srdb is None:
            data = datadir + "/data/SRSo.csv"
            srdb = srdb_from_csv(data)

        if geoLeaf is None:
            geoLeaf = genGeoLeaf()

        if dynamic_bins is None:
            self.dynamic = False
        else:
            self.dynamic = True

        self.xydb = xydb
        self.srdb = srdb
        self.geoLeaf = geoLeaf
        self.bins = dynamic_bins
        self.discretisation_level = discretisation_level
        self.twist = twist
        self.fit_leaves()

    def fit_leaves(self):
        leaves = {}
        xy = self.xydb
        sr = self.srdb
        for k in list(xy.keys()):
            leaves[k] = []
            for i in range(len(xy[k])):
                if self.dynamic:
                    xysr = {
                        age: (
                            list(xy[k][i][age]["x"]),
                            list(xy[k][i][age]["y"]),
                            list(sr[k]["s"]),
                            list(sr[k]["r"]),
                        )
                        for age in xy[k][i]
                    }
                else:
                    xysr = (xy[k][i][0], xy[k][i][1], sr[k][0], sr[k][1])
                leaves[k].append(xysr)
        leaves, discard = fitting.fit_leaves(
            leaves, self.discretisation_level, self.dynamic
        )

        self.leaves = leaves

    def get_age_index(self, age=None):
        age_index = age
        if age is not None:
            binf = (
                numpy.searchsorted(self.bins, age) - 1
            )  # binf is position of bin below age (starting at zero)
            if binf < 0:
                binf = 0  # age below first value are in firsts interval
            if binf >= (len(self.bins) - 1):  # binf is the last value
                age_index = binf - 1
            else:
                age_index = binf
        return age_index

    def get_leaf_key(self, lindex, lseed, age=None):
        # to do return one default leaf even if key error occur ?
        key, index = leaf_keys(lindex, lseed, self.srdb)
        age_index = self.get_age_index(age)
        # age_index = '(%s, %s]'%(str(self.bins[age_index-1]), str(self.bins[age_index]))
        return key, index, age_index

    def get_leaf(self, leaf_key):
        key, index, age_index = leaf_key
        if age_index is None:
            leaf = self.leaves[key][index]
        else:
            leaf = self.leaves[key][index][age_index]
        # deep copy is required as fitting.mesh alter shape
        leaf = deepcopy(leaf)
        if isinstance(leaf, dict):
            leaf = leaf["x"], leaf["y"], leaf["s"], leaf["r"]
        return leaf

    def get_sr(self, leaf_key):
        key, index, age_index = leaf_key
        sr = self.srdb[key]
        if isinstance(sr, dict):
            sr = sr["s"], sr["r"]
        return sr

    def blade_elt_area(self, leaf_key, Lshape, Lwshape, sr_base, sr_top):
        """surface of a blade element, positioned with two relative curvilinear abscissa"""

        S = 0
        sr_base = min([1, max([0, sr_base])])
        sr_top = min([1, max([sr_base, sr_top])])
        if leaf_key is not None:
            s, r = self.get_sr(leaf_key)
            sre = [sr for sr in zip(s, r) if (sr_base < sr[0] < sr_top)]
            if len(sre) > 0:
                se, re = list(zip(*sre))
                snew = [sr_base] + list(se) + [sr_top]
                rnew = (
                    [numpy.interp(sr_base, s, r)]
                    + list(re)
                    + [numpy.interp(sr_top, s, r)]
                )
            else:
                snew = [sr_base, sr_top]
                rnew = [numpy.interp(sr_base, s, r), numpy.interp(sr_top, s, r)]

            S = simpson(rnew, x=snew) * Lshape * Lwshape

        return S

    def blade_elt_width(self, leaf_key, Lshape, Lwshape, sr_base, sr_top):
        """width of a blade element, positioned with two relative curvilinear abscissa"""

        w = 0
        sr_base = min([1, max([0, sr_base])])
        sr_top = min([1, max([sr_base, sr_top])])
        if leaf_key is not None:
            s, r = self.get_sr(leaf_key)
            sre = [sr for sr in zip(s, r) if (sr_base < sr[0] < sr_top)]
            if len(sre) > 0:
                se, re = list(zip(*sre))
                rnew = (
                    [numpy.interp(sr_base, s, r)]
                    + list(re)
                    + [numpy.interp(sr_top, s, r)]
                )
            else:
                rnew = [numpy.interp(sr_base, s, r), numpy.interp(sr_top, s, r)]

            w = numpy.mean(rnew) * Lwshape

        return w

    def mesh(
        self, leaf_key, L_shape, Lw_shape, length, s_base, s_top, incline=1, flipx=False, min_area=1e-6
    ):
        """Compute mesh for a leaf element.
        - shape is a x,y,s,r tuple describing leaf shape
        - L_shape is the length of the scaled shape
        - Lw_shape is the width of the scaled shape
        - length is the total visible length to be meshed
        - s_base and s_top are relative proportion (on length) of the element to represent
        """
        shape = self.get_leaf(leaf_key)

        shape = incline_leaf(shape, incline)
        if flipx:
            shape = (-shape[0],) + shape[1:]  # to position leaves along tiller emitted
        leaf_mesh = fitting.mesh4(
            shape, L_shape, length, s_base, s_top, Lw_shape, twist=self.twist,
        min_area=min_area)
        if leaf_mesh:
            pts, ind = leaf_mesh
            if len(ind) < 1:
                # raise  AdelError('ERROR less than 1 triangles') # mesh4 filters triangle < 1e-6
                mesh = None
            else:
                mesh = fitting.plantgl_shape(pts, ind)
        else:
            if length > 0:
                print("ERROR No mesh", s_base, s_top, length)
                pass
            mesh = None

        return mesh
    def form_factor(self):
        """
        return form factor for each key in sr_db
        """
        try:
            return {
                k: simpson(self.srdb[k]["r"], x=self.srdb[k]["s"]) for k in self.srdb
            }
        except TypeError:
            return {k: simpson(self.srdb[k][1], x=self.srdb[k][0]) for k in self.srdb}

    def midrib(self, blade, resample=False):
        """Compute visible midrib x,y coordinates  and vertical distance to insertion point due to rolling from a blade node"""
        x = y = dy = None
        if blade.visible_length > 0.01:
            if blade.shape_key is not None:
                if self.dynamic:
                    inclin = 1  # inclination is encoded in db
                else:
                    inclin = blade.inclination

                shape = self.get_leaf(blade.shape_key)
                shape = incline_leaf(shape, inclin)

                x, y, s, r = fitting.leaf_element(
                    shape,
                    blade.shape_mature_length,
                    blade.visible_length,
                    0,
                    1,
                    blade.shape_max_width,
                )

                if resample:

                    def _resample(xx, yy):
                        ss = curvilinear_abscisse(xx, yy)
                        snew = numpy.linspace(0, max(ss), 20)
                        return numpy.interp(snew, ss, xx), numpy.interp(snew, ss, yy)

                    x, y = _resample(x, y)

                dy = blade.rolled_length

        return x, y, dy


def leaves_node(
    xydb=None,
    srdb=None,
    geoLeaf=None,
    dynamic_bins=None,
    discretisation_level=9,
    twist=0,
):
    return Leaves(**locals())
