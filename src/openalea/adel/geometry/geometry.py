from openalea.adel import fitting, symbol
import openalea.plantgl.all as pgl


def leaf_to_mesh(leaf, lmax, l, rmax, antisens):
    pts, ind = fitting.mesh3(leaf, lmax, l, rmax, antisens)
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh


def leaf_to_mesh_new(leaf, lmax, l, rmax, twist=False, nb_twist=1.0, nb_waves=8):
    pts, ind = fitting._mesh(
        leaf,
        lmax,
        l,
        rmax,
        functor=fitting.leaf_to_mesh_new,
        twist=twist,
        nb_twist=nb_twist,
        nb_waves=nb_waves,
    )
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh


def leaf_to_mesh_cicloid_twist(leaf, lmax, l, rmax, twist=True, nb_twist=2, nb_waves=8):
    pts, ind = fitting.mesh3(leaf, lmax, l, rmax)
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh


def leaf_element(leaf, lmax, l, s_base, s_top, rmax):
    pts, ind = fitting.mesh4(leaf, lmax, l, s_base, s_top, rmax)
    mesh = fitting.plantgl_shape(pts, ind)
    return mesh


def symbols(database, seed=None, relative_angle=True):
    return symbol.build_symbols(database, seed, relative_angle=relative_angle)


def LeafElement(sym, leaf_rank, length, final_length, radius_max, incB, index):
    leaf = sym.get("LeafElement")
    element = leaf(
        1, final_length, length, radius_max, 0, 1, leaf_rank, 1, incB, index - 1
    )
    return element["geometry"], element["label"]


def mesh2shapes(scene):
    """Convert all the meshes containing colors into independent shapes."""
    if isinstance(scene, pgl.Scene):
        scene = scene
    else:
        # Case caribu scene
        scene = scene.scene

    new_scene = pgl.Scene()
    for shape in scene:
        g = shape.geometry
        if isinstance(g, pgl.TriangleSet) and len(g.colorList) > 0:
            # convert the mesh into shapes with different colors
            pts = g.pointList
            indx = g.indexList
            colors = g.colorList
            for i, ind in enumerate(indx):
                new_geom = pgl.TriangleSet([], [])
                new_geom.indexList.append(pgl.Index3(0, 1, 2))
                for k in ind:
                    new_geom.pointList.append(pts[k])
                _shape = pgl.Shape(new_geom, pgl.Material(pgl.Color3(colors[i])))
                new_scene.add(_shape)
        else:
            new_scene.add(shape)
    return (new_scene,)
