import numpy
from scipy.interpolate import splev, splprep
import openalea.adel.fitting as fitting
from openalea.adel.data_samples import leaves_db

# from openalea.plantgl.all import Viewer
import openalea.plantgl.all as pgl

# nosetetsts fails importing pylab for some mysteriuous reason (backend ?)
with_pylab = True
try:
    from pylab import plot, clf
except:
    with_pylab = False


leaves = leaves_db()
rank = list(leaves.keys())[0]
leaf = leaves[rank][0]

#
# build a database of function
# 1. fit leaves
# 2. write a function
# 3. draw the result

# 3 bis. write a class
# 4 parse the string and build the scene
# 5. create the can file
translation = pgl.Vector3(0, 0, 0)
zt = 6
yt = 13


def test1(leaf=leaf, scene=None):
    """
    Fit the leaf, create a mesh and draw it with plantGL.
    """
    global translation
    if scene is None:
        scene = pgl.Scene()
        # Viewer.display(scene)

    x, y, s, r = leaf
    spline_leaf, leaf_surface = fitting.fit_leaf(x, y, s, r)
    mesh = fitting.discretize(spline_leaf, 30, 7, 1)

    scene += pgl.Translated(translation, mesh)
    # Viewer.update()


# test1(leaf)
def test2(leaves=leaves):
    """
    Visualize all the leaves in the database in the same scene.
    """
    global translation, yt, zt
    translation = pgl.Vector3(0, 0, 0)
    scene = pgl.Scene()
    # Viewer.display(scene)
    for k in leaves:
        print("Rank number: ", k)
        index = 0
        for leaf in leaves[k]:
            try:
                test1(leaf, scene)
            except:
                print("problem with leaf %d in rank %s" % (index, k))
            index += 1
            translation.z += zt
            # raw_input('Enter something for next leaf...')
        translation.y = -60 + yt * int(k)
        translation.z = 0


def test3(leaf=leaf, scene=None):
    """
    Obtain a leaf at different time.
    """
    global zt

    if scene is None:
        scene = pgl.Scene()

    x, y, s, r = leaf
    spline_leaf, leaf_surface = fitting.fit_leaf(x, y, s, r)
    translation = pgl.Vector3()
    n = 7

    for i in range(1, n + 1):
        mesh = fitting.partial_leaf(spline_leaf, 30, n, i, 1)
        scene += pgl.Translated(translation, mesh)
        translation.z += zt

    # Viewer.display(scene)


# def test4(leaves=leaves,rank=rank):
#     # try with rank=7
#     index = 0
#     if rank in leaves:
#         for leaf in leaves[rank]:
#             test3(leaf)
#             raw_input('rank %s index %d'%(rank,index))
#             index += 1


def test5(leaf=leaf):
    """
    Fit the leaf, create a mesh and draw it with plantGL.
    """
    x, y, s, r = leaf
    # spline_leaf, leaf_surface = fitting.fit_leaf(x, y, s, r)
    tckp, u = splprep([x, y], k=5)
    xp, yp = splev(numpy.linspace(0, 1, 6), tckp)
    # xn, yn, rn = splev(linspace(0,1,100),spline_leaf)
    sn = fitting.curvilinear_abscisse(x, y)
    l = sn.sum()
    if with_pylab:
        clf()
        plot(x, y, ".")
        plot(xp, yp)


# def test6(leaves=leaves,rank=rank):
#     # try with rank=7
#     index = 0
#     if rank in leaves:
#         for leaf in leaves[rank]:
#             test5(leaf)
#             q=raw_input('rank %s index %d'%(rank,index))
#             if q=='q':
#                 break
#             index += 1


def test7(leaf=leaf, scene=None):
    if scene is None:
        scene = pgl.Scene()
        # Viewer.display(scene)

    x, y, s, r = leaf
    spline_leaf, leaf_surface = fitting.fit_leaf(x, y, s, r)
    pts, ind = fitting.mesh(spline_leaf, 30, 7, 7, 1)
    # fitting.write_smf('leaf_full.smf', pts, ind)
    # Viewer.display(fitting.plantgl_shape(pts, ind))


def test8(leaf=leaf, scene=None):
    global translation, zt
    if scene is None:
        scene = pgl.Scene()
        # Viewer.display(scene)

    x, y, s, r = leaf
    leaf_new, leaf_surface = fitting.fit2(x, y, s, r)

    pts, ind = fitting.mesh2(leaf_new, 7, 7, 1)
    # pts2, ind2 = fitting.qslim(13, pts, ind)
    # mesh = fitting.plantgl_shape(pts2, ind2)

    # sc=pgl.SurfComputer(pgl.Discretizer())
    # mesh.apply(sc)
    # scale_z = leaf_surface*7 / (sc.surface)
    # mesh_final = mesh.transform(pgl.Scaling((1,1,scale_z)))
    # mesh_final = mesh
    scene += pgl.Translated(translation, fitting.plantgl_shape(pts, ind))
    # scene += pgl.Translated(translation+(0,yt/3.,0), mesh_final)

    # Viewer.update()


def test81(leaf=leaf, scene=None):
    global translation, yt, zt
    if scene is None:
        scene = pgl.Scene()
        # Viewer.display(scene)

    mesh = fitting.leaf_shape2(leaf, 10, 7, 7, 1)

    scene += pgl.Translated(translation + (0, yt / 3.0, 0), mesh)

    # Viewer.update()


# def test9(leaves=leaves):
#     """
#     Visualize all the leaves in the database in the same scene.
#     """
#     global translation, yt, zt
#     translation = pgl.Vector3(0,0,0)
#     scene= pgl.Scene()
#     Viewer.display(scene)
#     for k in leaves:
#         print "Rank number: ", k
#         index = 0
#         for leaf in leaves[k]:
#             try:
#                 test8(leaf,scene)
#                 test81(leaf, scene)
#             except:
#                 print "problem with leaf %d in rank %s"%(index,k)
#             index += 1
#             translation.y += yt
#             #raw_input('Enter something for next leaf...')
#         translation.z = -60 + zt*int(k)
#         translation.y = 0
#
#
# def test10():
#     pass
# sc=pgl.SurfComputer(pgl.Discretizer())
# m=pgl_shape mesh2
# m2=pgl_shape qslim
# m.apply(sc)
# sc.surface
# m2.apply(sc)
# sc.surface
