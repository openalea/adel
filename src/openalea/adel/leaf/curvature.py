import numpy as np

# from openalea.plantgl.all import *
from .curve_discretizer import curve_discretizer
# from openalea.adel.fitting import curvilinear_abscisse


def curvilinear_abscisse(x, y):
    s = np.zeros(len(x))
    s[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    return s.cumsum()


def curvature_xys(x, y, s):
    ds = np.diff(s)
    dx, dy = np.diff(x), np.diff(y)
    # dx /= ds
    # dy /= ds
    theta = np.arctan2(dy, dx)

    dtheta = np.diff(theta) / ds[1:]
    return (x[0], y[0]), theta[0], s, dtheta


def curvature(crv, n=100):
    """Compute the curvature of a 2D curve.

    Return the first point, the first tangent, the curvilinear abscissa and
    the curvature.
    """

    x, y = curve_discretizer(crv, n)
    s = curvilinear_abscisse(x, y)
    L = s.max()
    s /= L
    x /= L
    y /= L
    ds = np.diff(s)
    dx, dy = np.diff(x), np.diff(y)
    # dx /= ds
    # dy /= ds
    theta = np.arctan2(dy, dx)

    dtheta = np.diff(theta) / ds[1:]
    return (x[0], y[0]), theta[0], s, dtheta


def curvature2xy(p0, angle, s, dtheta):
    x0, y0 = p0
    ds = np.diff(s)
    theta = angle + np.cumsum([0] + list(dtheta * ds[1:]))

    dx = ds * np.cos(theta)
    dy = ds * np.sin(theta)

    x = np.cumsum([0] + list(dx)) + x0
    y = np.cumsum([0] + list(dy)) + y0

    return x, y


def interpolate_curvature(curvatures, times, kind="cubic"):
    """Interpolate the curvatures.

    A curvature is a parametrisation `s`, d(angle)/ ds and a parameter between [0,1].
    Return a surface f(s,t).
    """

    from scipy import interpolate
    from scipy.ndimage import measurements

    if len(curvatures) == 3 and not isinstance(curvatures[2], tuple):
        curvatures = [curvatures]

    curv_abs = [c[0] for c in curvatures]
    curves = [c[1] for c in curvatures]
    params = [c[2] for c in curvatures]

    n = len(params)
    if n == 1:
        curv_abs *= 2
        curves *= 2
        curves[0] = np.zeros(len(curves[0]))
        params = [0.0, 1.0]
    # elif False:
    #     if params[0] != 0.0:
    #         params.insert(0, 0.0)
    #         curves.insert(0, curves[0])
    #         curv_abs.insert(0, curv_abs[0])
    #     if params[-1] != 0:
    #         params.append(1.0)
    #         curves.append(curves[-1])
    #         curv_abs.append(curv_abs[-1])

    # compute a common parametrisation s
    # We conserve the last parametrisation because the result is very sensitive
    if True:
        min_s = min(np.diff(s).min() for s in curv_abs)
        s_new = np.unique(np.array(curv_abs).flatten())
        ds = np.diff(s_new)
        k = np.cumsum(ds >= min_s)
        labels = list(range(k.min(), k.max() + 1))
        s = np.zeros(len(labels) + 1)
        s[1:] = measurements.mean(s_new[1:], k, labels)
        s[-1] = 1.0
        s = s_new
    # else:
    #     s = np.array(curv_abs[-1])

    # renormalise all the curves
    curves = [
        np.interp(s[1:-1], old_s[1:-1], old_crv)
        for old_s, old_crv in zip(curv_abs, curves)
    ]

    # interpolate correctly the curvatures
    x = s[1:-1]
    y = np.array(params)
    z = np.array(curves)

    # f = interpolate.interp2d(y, x, z, kind=kind)
    f = interpolate.RectBivariateSpline(x, y, z.T, kx=1, ky=1)
    return s, f(x, times)


def curvatures2xy(p0, angle, s, curvatures, index):
    return curvature2xy(p0, angle, s, curvatures[:, index])
