import os
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.integrate import simpson, trapezoid


import openalea.plantgl.all as pgl
from .simplification import cost


##### DEBUG
debug = False


def curvilinear_abscisse(x, y):
    s = np.zeros(len(x))
    s[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    return s.cumsum()


def fit_leaf(x, y, s, r):
    """
    x, y: 2d points
    s: curvilinear abscissa
    r: leaf radius

    Algo:
        #. fit x, y
            #. discretize the splines (high: 100)
            #. compute curvilinear abscissa
            #. normalize :
                + compute length and divide x, y by length

        #. fit r, s
            #. discretize r, s (high 500)
            #. normalize :
                + compute max r, s and divide r, s
            #. compute leaf surface
            #. r-t1-> found r(s): interp(s_t1,s_t2, r_t2)

        #. fit x(t1), y(t1), r(t1)
        #. return (x, y, r) and surface value

    """
    global debug
    # spline parameters
    smooth = 3.0  # smoothness parameter
    k = 3  # spline order
    nest = -1  # estimate of number of knots needed (-1 = maximal)

    if debug:
        tckp, u = splprep([x, y])
    else:
        tckp, u = splprep([x, y], s=smooth, k=k, nest=nest)

    xnew, ynew = splev(np.linspace(0, 1, 100), tckp)
    snew = curvilinear_abscisse(xnew, ynew)
    # 1.4
    snew = snew / snew.max()

    # 2.1
    if debug:
        tckp2, v = splprep([s, r])
    else:
        tckp2, v = splprep([s, r], s=smooth, k=k, nest=nest)
    snew2, rnew2 = splev(np.linspace(0, 1, 500), tckp2)

    snew2 = snew2 / snew2.max()
    rnew2 = rnew2 / rnew2.max()

    # 2.4 leaf surface (integral r ds)
    leaf_surface = simpson(rnew2, x=snew2)

    # 2.5 Compute rnew
    rnew = np.interp(snew, snew2, rnew2)
    if debug:
        # plot( s/s.max(), r/r.max())
        # plot( snew, rnew )
        pass
    # radius always positive
    # rnew = where(rnew>0., rnew, zeros(len(rnew)))

    # 3.
    if debug:
        tckp_final, t = splprep([xnew, ynew, rnew])
        # xf, yf, rf = splev(linspace(0,1,15),tckp_final)
        # plot(x/x.max(), y/x.max())
        # plot (xf/xf.max(), yf/xf.max())
        # plot(xf/xf.max(), rf/xf.max())
    else:
        tckp_final, t = splprep([xnew, ynew, rnew], s=smooth, k=k, nest=nest)
    # xf, yf, rf = splev(linspace(0,1,15),tckp_final)
    return tckp_final, leaf_surface


def fit2(x, y, s, r):
    """
    x, y: 2d points
    s: curvilinear abscissa
    r: leaf radius

    Algo:
        #. fit x, y
            #. discretize the splines (high: 100)
            #. compute curvilinear abscissa
            #. normalize :
                + compute length and divide x, y by length

        #. fit r, s
            #. discretize r, s (high 500)
            #. normalize :
                + compute max r, s and divide r, s
            #. compute leaf surface
            #. r_t1-> found r(s): interp(s_t1,s_t2, r_t2)

        #. fit x(t1), y(t1), r(t1)
        #. return (x, y, r) and surface value

    """
    global debug
    # spline parameters
    smooth = len(x) - np.sqrt(2 * len(x))  # smoothness parameter
    smooth = 0.01

    k = 3  # spline order
    nest = -1  # estimate of number of knots needed (-1 = maximal)

    # if debug:
    #    tckp,u = splprep([x,y])
    # else:
    #    tckp,u = splprep([x,y],s=smooth,k=k,nest=nest)

    try:
        tckp, u = splprep([x, y], s=smooth, k=k, nest=nest)
    except:
        try:
            tckp, u = splprep([x, y])
        except:
            tckp, u = splprep([x, y], k=1)

    xnew, ynew = splev(np.linspace(0, 1, 100), tckp)
    snew = curvilinear_abscisse(xnew, ynew)
    # 1.4
    length = snew.max()
    snew /= length
    xnew /= length
    ynew /= length

    # 2.1
    try:
        tckp2, v = splprep([s, r], s=smooth / 100.0, k=k, nest=nest)
    except:
        tckp2, v = splprep([s, r])

    snew2, rnew2 = splev(np.linspace(0, 1, 200), tckp2)

    snew2 = snew2 / snew2.max()
    rnew2 = np.maximum(rnew2, 0) / rnew2.max()

    # 2.4 leaf surface (integral r ds)
    leaf_surface = simpson(rnew2, x=snew2)

    # 2.5 Compute rnew
    rnew = np.interp(snew, snew2, rnew2)
    return (xnew, ynew, snew, rnew), leaf_surface


def discretize(leaf_spline, nb_polygones, length_max, radius_max):
    """
    Compute a mesh from a leaf represented by a 3d spline containing
    x, y and radius.

    Final radius have to be null.
    """
    return partial_leaf(leaf_spline, nb_polygones, length_max, length_max, radius_max)


def partial_leaf(leaf_spline, nb_polygones, length_max, length, radius_max):
    """
    Compute a mesh from a leaf represented by a 3d spline containing
    x, y and radius.
    Final radius have to be null.
    The leaf is obtained by evaluating the central curve on a domain
    [0, l/lmax] and the radius on a domain [1-l/lmax, 1].

    """
    if length <= 0.0:
        return
    if length > length_max:
        length = length_max

    param = float(length) / length_max

    xf, yf, rnull = splev(np.linspace(0, param, nb_polygones), leaf_spline)
    xnull, ynull, rf = splev(np.linspace(1 - param, 1, nb_polygones), leaf_spline)
    rf = np.where(rf > 0.0, rf, np.zeros(len(rf)))
    rf[-1] = 0.0

    xf *= length_max
    yf *= length_max
    rf *= radius_max

    # build a mesh
    n = len(xf)
    points = list(zip(xf, yf, -rf / 2.0))
    points.extend(list(zip(xf, yf, rf / 2.0)))

    ind = np.array(range(n - 2))
    indices = list(zip(ind, ind + n, ind + (n + 1)))
    indices.extend(list(zip(ind, ind + (n + 1), ind + 1)))
    # add only one triangle at the end !!
    indices.append((n - 2, 2 * n - 2, n - 1))

    return plantgl_shape(points, indices)


def mesh(leaf_spline, nb_polygones, length_max, length, radius_max):
    if length <= 0.0:
        return
    if length > length_max:
        length = length_max

    param = float(length) / length_max

    xf, yf, rnull = splev(np.linspace(0, param, nb_polygones), leaf_spline)
    xnull, ynull, rf = splev(np.linspace(1 - param, 1, nb_polygones), leaf_spline)
    rf = np.where(rf > 0.0, rf, np.zeros(len(rf)))
    rf[-1] = 0.0
    xf *= length_max
    yf *= length_max
    rf *= radius_max
    return leaf_to_mesh(xf, yf, rf)


def mesh3(leaf, length_max, length, radius_max, antisens=True):
    return _mesh(leaf, length_max, length, radius_max, antisens=antisens)


def leaf_to_mesh_new(x, y, r, twist=True, nb_twist=1.0, nb_waves=8, **kwds):
    # twist = True
    # nb_twist = 1.

    # test sin
    # 7.64 is the length of a sin arc from 0 to 2pi
    s = curvilinear_abscisse(x, y)
    L = s.max()
    s /= L

    nb_oscillation = nb_waves
    dt = list(np.arctan2(np.diff(y), np.diff(x)))
    dt.append(0.0)
    dt = np.array(dt)
    if twist:
        angle = 2 * np.pi * nb_twist * s
    else:
        angle = np.pi * nb_oscillation * s
    if not twist:
        waves1_x = 0
        waves1_y = np.sin(angle) * r / 2.0 * s
        waves2_x = 0
        waves2_y = np.sin(angle + np.pi) * r / 2.0 * s
        waves_z = 1
    else:
        scalar = r / 2.0
        waves1_y = -np.sin(angle)
        waves1_x = 0 * scalar
        waves2_x = 0 * scalar
        waves2_y = np.sin(angle)
        waves_z = np.cos(angle)

    n = len(x)
    points = list(zip(x + waves1_x, -r / 2.0 * waves_z, y + waves1_y))
    points.extend(list(zip(x, np.zeros(n), y)))
    points.extend(list(zip(x + waves2_x, r / 2.0 * waves_z, y + waves2_y)))

    ind = np.array(range(n - 2))
    indices = list(zip(ind, ind + n, ind + (n + 1)))
    indices.extend(list(zip(ind, ind + (n + 1), ind + 1)))
    indices.extend(list(zip(ind + n, ind + 2 * n, ind + 2 * n + 1)))
    indices.extend(list(zip(ind + n, ind + 2 * n + 1, ind + n + 1)))

    # add only one triangle at the end !!
    if r[-1] < 0.001:
        indices.append((n - 2, 3 * n - 2, 2 * n - 1))
    else:
        indices.append((n - 2, 2 * n - 2, 2 * n - 1))
        indices.append((n - 2, 2 * n - 1, n - 1))
        indices.append((2 * n - 2, 3 * n - 2, 3 * n - 1))
        indices.append((2 * n - 2, 3 * n - 1, 2 * n - 1))

    return points, indices


def leaf_to_mesh(x, y, r, twist_start=0, twist_end=0, **kwds):
    n = len(x)
    theta = np.linspace(np.radians(twist_start), np.radians(twist_end), n)
    points = list(
        zip(x, -r / 2.0 * abs(np.cos(theta)), y + abs(np.sin(theta)) * r / 2.0)
    )
    points.extend(
        list(zip(x, r / 2.0 * abs(np.cos(theta)), y - abs(np.sin(theta)) * r / 2.0))
    )

    ind = np.array(range(n - 2)) if n > 2 else np.array([0])
    indices = list(zip(ind, ind + n, ind + (n + 1)))
    indices.extend(list(zip(ind, ind + (n + 1), ind + 1)))

    # add only one triangle at the end !!
    if n < 2:
        assert len(points) == 4
        if r[-1] < 0.001:
            indices = indices[0:1]
    elif r[-1] < 0.001:
        indices.append((n - 2, 2 * n - 2, 2 * n - 1))
    else:
        indices.append((n - 2, 2 * n - 2, 2 * n - 1))
        indices.append((n - 2, 2 * n - 1, n - 1))

    return points, indices


def _mesh(
    leaf, length_max, length, radius_max, antisens=True, functor=leaf_to_mesh, **kwds
):
    from openalea.adel.leaf.curvature import curvature_xys, curvature2xy

    if length <= 0.0:
        return
    if length > length_max:
        length = length_max

    x, y, s, r = leaf
    param = float(length) / length_max

    n = len(x)

    sample_sr = np.linspace(1.0 - param, 1.0, num=n)
    if antisens:
        sample_xy = np.linspace(0, param, num=n)
    else:
        sample_xy = sample_sr

    xn = np.interp(sample_xy, s, x) * length_max
    yn = np.interp(sample_xy, s, y) * length_max
    if not antisens:
        p, theta, _s, dtheta = curvature_xys(x, y, s)
        pn, thetan, sn, dthetan = curvature_xys(xn, yn, sample_xy * length_max)
        xn, yn = curvature2xy((0, 0), theta, sn, dthetan)
    rn = np.interp(sample_sr, s, r) * radius_max

    return functor(xn, yn, rn, **kwds)


def mesh2(leaf, length_max, length, radius_max):
    return mesh4(leaf, length_max, length, 0.0, 1.0, radius_max)


def leaf_element(leaf, length_max, length, s_base, s_top, radius_max):
    def insert_values(a, values):
        l = a.tolist()
        l.extend(values)
        return np.unique(l)

    s_base = min(s_base, s_top, 1.0)
    s_top = max(s_base, s_top, 0.0)

    if length <= 0.0:
        return
    if length > length_max:
        length = length_max

    x, y, s, r = leaf
    # just scale leaf case
    if length == length_max and s_base == 0 and s_top == 1:
        x *= length_max
        y *= length_max
        s *= length_max
        r *= radius_max
        return x, y, s, r

    # 1. compute s_xy and s_r for length vs length_max

    # force the leaf width to zero at the top
    r[-1] = 0

    param = float(length) / length_max
    n = len(s)
    # add param in s_xy
    sp = insert_values(s, [1 - param, param])
    s_xy = np.compress(sp <= param, sp)
    s_r = np.compress(sp >= (1 - param), sp)
    s_xy = np.union1d(s_xy, s_r - s_r[0])
    s_r = s_xy + s_r[0]
    # add s_base and s_top in s

    s_base *= param
    s_top *= param
    s_valid = insert_values(s_xy, [s_base, s_top])
    s_valid = np.compress(s_valid <= s_top, s_valid)
    s_valid = np.compress(s_valid >= s_base, s_valid)

    # delete small intervals COMIT from  Here !
    eps = (s_top - s_base) / (len(s) * 2)
    ds = s_valid[1:] - s_valid[:-1]
    error = ds >= eps
    s_val = np.compress(error, s_valid[1:])
    s_val = insert_values(s_val, [s_valid[0], s_valid[-1]])

    # find param and 1-param in s...
    # build xf, yf, rf with the same nb of points

    s_xyf = s_val
    s_rf = s_val + s_r[0]

    xf = np.interp(s_xyf, s, x)
    yf = np.interp(s_xyf, s, y)
    rf = np.interp(s_rf, s, r)

    cond = rf <= 0
    if cond.any():
        # delete points with negative radius
        index = cond.searchsorted(True)
        xf = xf[: index + 1]
        yf = yf[: index + 1]
        rf = rf[: index + 1]
        rf[-1] = 0.0

    xf *= length_max
    yf *= length_max
    rf *= radius_max

    return xf, yf, s_val, rf


def mesh4(leaf, length_max, length, s_base, s_top, radius_max, twist=0, volume=0.1,min_area=1e-6):
    xf, yf, s_val, rf = leaf_element(
        leaf, length_max, length, s_base, s_top, radius_max
    )

    if len(xf) < 2:
        # All the radius are negative or null.
        # Degenerated element.
        return [], []

    pts, ind = leaf_to_mesh(
        xf, yf, rf, twist_start=twist * min(s_val), twist_end=twist * max(s_val)
    )

    def _surf(ind, pts):
        from openalea.plantgl.all import norm, cross, Vector3

        A, B, C = [Vector3(pts[i]) for i in ind]
        return norm(cross(B - A, C - A)) / 2.0

    ind = [id for id in ind if _surf(id, pts) > min_area]
    return pts, ind


def write_smf(filename, points, indices):
    """
    Write a smf file for QSlim software.
    See http://people.scs.fsu.edu/~burkardt/data/smf/smf.html
    """

    header = """#$SMF 1.0
#$vertices %d
#faces %d
#

""" % (len(points), len(indices))

    pts_str = "\n".join(["v %f %f %f" % (pt[0], pt[1], pt[2]) for pt in points])
    ind_str = "\n".join(
        ["f %d %d %d" % (ind[0] + 1, ind[1] + 1, ind[2] + 1) for ind in indices]
    )

    f = open(filename, "w")
    f.write(header)
    f.write(pts_str)
    f.write("\n")
    f.write(ind_str)

    f.close()


def read_smf(filename):
    """
    Read a smf file containing a mesh.
    Returns the point set and the face set.
    """
    f = open(filename)
    points = []
    indices = []

    for l in f:
        l.strip()
        if not l:
            continue
        if l[0] == "v":
            pts = l.split()[1:]
            points.append([float(pt) for pt in pts])
        elif l[0] in ["f", "t"]:
            ind = l.split()[1:]
            indices.append([int(i) - 1 for i in ind])
        else:
            continue

    f.close()

    return points, indices


def plantgl_shape(points, indices):
    indices = [list(map(int, index)) for index in indices]
    return pgl.TriangleSet(points, indices, normalPerVertex=False)


def qslim(nb_triangles, points, indexes):
    """
    Decimate a mesh using the qslim software.
    Return the new mesh with fewer triangles.
    Try to respect angles and end points.
    """

    file_in = os.tempnam() + ".smf"
    file_out = os.tempnam() + ".smf"

    write_smf(file_in, points, indexes)

    pts, ind = read_smf(file_out)

    os.remove(file_in)
    os.remove(file_out)

    return pts, ind


def leaf_shape(leaf, nb_triangles, length_max, length, radius_max):
    if length <= 0:
        return
    x, y, s, r = leaf
    leaf_new, leaf_surface = fit2(x, y, s, r)

    pts, ind = mesh2(leaf_new, length_max, length, radius_max)
    if len(pts) <= 1.5 * nb_triangles:
        pts2, ind2 = pts, ind
    else:
        pts2, ind2 = qslim(nb_triangles, pts, ind)
    # pts2, ind2 = pts, ind
    mesh = plantgl_shape(pts2, ind2)

    sc = pgl.SurfComputer(pgl.Discretizer())
    mesh.apply(sc)
    scale_radius = leaf_surface * length_max * radius_max / (sc.surface)
    _ = mesh.transform(pgl.Scaling((1, scale_radius, 1)))
    mesh_final = mesh
    return mesh_final


def fit3(x, y, s, r, nb_points):
    leaf, leaf_surface = fit2(x, y, s, r)  # use here simpson as leaf is a smooth shape
    xn, yn, sn, rn = simplify(leaf, nb_points)
    return xn, yn, sn, rn


def simplify(leaf, nb_points, scale_radius=True):
    """ " Simplify a 2d polyline up to nb_points
    Parameters
    ----------
    leaf : a x, y, s, r tuple of iterables representing the 2d polyline to be simplified
    nb_points : the number of points along the polyline after simplification
    scale_radius : (bool, default=True)  should radius be scaled after simplification to ensure
        polygonial_area(simplified_leaf) = smooth_area(input_leaf) ?

    Returns
    -------
    a x, y, s, r tuple of iterables for the simplified 2d polyline
    """
    xn, yn, sn, rn = leaf

    points = [pgl.Vector3(*pt) for pt in zip(xn, rn, yn)]
    pts = [_f for _f in cost(points, nb_points) if _f]
    coords = ((pt.x, pt.y, pt.z) for pt in pts)
    x, r, y = list(map(np.array, zip(*coords)))
    s = curvilinear_abscisse(x, y)
    # keep smax similar to sn
    adj = max(sn) / max(s)
    s *= adj
    x *= adj
    y *= adj

    if scale_radius:
        leaf_surface = simpson(rn, x=sn)  # use here simpson as leaf is a smooth shape
        new_surface = trapezoid(
            r, x=s
        )  # use here trapezoid as the new surface is a polygonal shape
        scale_r = leaf_surface / new_surface
        r *= scale_r

    return x, y, s, r


def leaf_shape2(leaf, nb_triangles, length_max, length, radius_max):
    if length <= 0:
        return

    x, y, s, r = leaf
    nb_points = nb_triangles / 2 + 2
    leaf_new = fit3(x, y, s, r, nb_points)

    pts, ind = mesh2(leaf_new, length_max, length, radius_max)
    mesh = plantgl_shape(pts, ind)

    return mesh


def _fit_element(el, nb_points):
    leaf = None
    if isinstance(el, dict):
        x, y, s, r = el["x"], el["y"], el["s"], el["r"]
    else:
        x, y, s, r = el
    try:
        leaf = fit3(x, y, s, r, nb_points)
        # force leaf tip to be s,r = 1,0
        leaf[2][-1] = 1.0
        leaf[3][-1] = 0.0
    except:
        pass
    return leaf


def fit_leaves(leaves, nb_points, dynamic=False):
    new_db = {}
    discarded = {}
    db = leaves

    for key in db:
        l = db[key]
        for i, el in enumerate(l):
            if not dynamic:
                leaf = _fit_element(el, nb_points)
            else:
                leaf = {age: _fit_element(v, nb_points) for age, v in el.items()}
                if any([x is None for x in list(leaf.values())]):
                    leaf = None
            if leaf is not None:
                new_db.setdefault(key, []).append(leaf)
            else:
                print(
                    (
                        "openalea.adel.fitting->fit_leaves: can't fit leaf shape index %s,"
                        " Rsub-index %d (python sub-index %d)=> leaf shape discarded"
                        % (key, i + 1, i)
                    )
                )
                discarded.setdefault(key, []).append(i + 1)

    return (
        new_db,
        discarded,
    )
