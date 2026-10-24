import numpy as np
cm_import_okay = True
try:
    from matplotlib import cm
    from matplotlib.colors import Normalize, LogNorm
except:
    cm_import_okay = False


""" A colormap class
    copied from  openalea.color.colormap because
    import openalea.color.colormap makes matplotlib
    with TkAgg frontend crash for some mysterious reason !!
"""


class ColorMap:
    """A RGB color map, between 2 colors defined in HSV code

    :Examples:

    >>> minh,maxh = minandmax([height(i) for i in s2])
    >>> colormap = ColorMap(minh,maxh)
    >>> s3 = [ Shape(i.geometry, Material
    >>>    (Color3(colormap(height(i))), 1), i.id)
    >>>    for i in s2]

    """

    def __init__(self, minval=0.0, maxval=1.0):
        self.minval = float(minval)
        self.maxval = float(maxval)

    def color(self, normedU):
        """

        :param normedU: todo

        """
        inter = 1 / 5.0
        winter = int(normedU / inter)
        a = (normedU % inter) / inter
        b = 1 - a

        if winter < 0:
            col = (self.coul2, self.coul2, self.coul1)
        elif winter == 0:
            col = (self.coul2, self.coul2 * b + self.coul1 * a, self.coul1)
        elif winter == 1:
            col = (self.coul2, self.coul1, self.coul1 * b + self.coul2 * a)
        elif winter == 2:
            col = (self.coul2 * b + self.coul1 * a, self.coul1, self.coul2)
        elif winter == 3:
            col = (self.coul1, self.coul1 * b + self.coul2 * a, self.coul2)
        elif winter > 3:
            col = (self.coul1, self.coul2, self.coul2)
        return int(col[0]), int(col[1]), int(col[2])

    def greycolor(self, normedU):
        """

        :param normedU: todo
        :returns: todo
        """
        return int(255 * normedU), int(255 * normedU), int(255 * normedU)

    def grey(self, u):
        """
        :param u:
        :returns: todo
        """
        return self.greycolor(self.normU(u))

    def normU(self, u):
        """
        :param u:
        :returns: todo
        """
        if self.minval == self.maxval:
            return 0.5
        return (u - self.minval) / (self.maxval - self.minval)

    def __call__(self, u, minval=0, maxval=1, coul1=80, coul2=20):
        self.coul1 = coul1
        self.coul2 = coul2
        self.minval = float(minval)
        self.maxval = float(maxval)
        return self.color(self.normU(u))


def get_cmap(cmap):
    _cmap = cmap
    if isinstance(cmap, str):
        _cmap = cm.get_cmap(cmap)
    return _cmap


def colormap(g, property_name, cmap="jet", lognorm=True):
    """Compute a color property based on a given property and a colormap."""
    prop = g.property(property_name)
    keys = list(prop.keys())
    values = np.array(list(prop.values()))
    # m, M = int(values.min()), int(values.max())

    if cm_import_okay:
        _cmap = get_cmap(cmap)
        norm = Normalize() if not lognorm else LogNorm()
        values = norm(values)
        # my_colorbar(values, _cmap, norm)

        colors = (_cmap(values)[:, 0:3]) * 255
        colors = np.array(colors, dtype=np.intp).tolist()
    else:
        cmap = ColorMap(values.min(), values.max())
        colors = cmap(values)

    g.properties()["color"] = dict(list(zip(keys, colors)))
    return g
