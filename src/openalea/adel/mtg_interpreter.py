"""Methods for mtg interpretation with turtle"""

from math import degrees, pi, cos, sin
import openalea.plantgl.all as pgl

# from openalea.mtg import MTG
from openalea.mtg.plantframe.turtle import TurtleFrame
# from openalea.mtg.algo import union
# from openalea.astk.plantgl_utils import addSets


def _is_iterable(x):
    try:
        x = iter(x)
    except TypeError:
        return False
    return True


def as_tuples(pgl_3List, offset=0):
    """return pgl list of 3 numbers kind (indes3, vector3) as a list of python tuples"""
    if not _is_iterable(offset):
        offset = [offset] * 3
    return [(i[0] + offset[0], i[1] + offset[1], i[2] + offset[2]) for i in pgl_3List]


def addSets(pglset1, pglset2, translate=(0, 0, 0)):
    """create a new TriangleSet by addition of two existing ones
    if translate is not None, pglset2 is translated with vector translate
    """
    points = as_tuples(pglset1.pointList) + as_tuples(
        pglset2.pointList, offset=translate
    )
    index = as_tuples(pglset1.indexList) + as_tuples(
        pglset2.indexList, offset=len(pglset1.pointList)
    )
    return pgl.TriangleSet(points, index)


# Meshing function for StemElements


def slim_cylinder(length, radius_base, radius_top):
    "Try to construct a cylinder with a low number of triangles."

    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2 * pi / 3.0, 4 * pi / 3.0
    r = rb
    p1 = (r * cos(a1), r * sin(a1), 0)
    p2 = (r * cos(a2), r * sin(a2), 0)
    p3 = (r * cos(a3), r * sin(a3), 0)
    r = rt
    q1 = (r * cos(a1 + pi), r * sin(a1 + pi), length)
    q2 = (r * cos(a2 + pi), r * sin(a2 + pi), length)
    q3 = (r * cos(a3 + pi), r * sin(a3 + pi), length)
    set = pgl.TriangleSet(
        [p1, p2, p3, q1, q2, q3],
        [
            (2, 1, 0),
            (3, 4, 5),
            (0, 5, 4),
            (0, 4, 2),
            (2, 4, 3),
            (3, 1, 2),
            (1, 3, 5),
            (5, 0, 1),
        ],
    )
    return set


def StemElement_mesh(length, diameter_base, diameter_top, classic=False):
    """Compute mesh for a stem element
    - classic indicates
    """
    if classic:
        solid = True
        diameter = diameter_base
        slices = 6  # 6 is the minimal number of slices for a correct computation of star (percentage error lower than 5)
        stem = pgl.Tapered(
            diameter_base / 2.0,
            diameter_top / 2.0,
            pgl.Cylinder(1.0, length, solid, slices),
        )
        tessel = pgl.Tesselator()
        stem.apply(tessel)
        mesh = tessel.triangulation
    else:
        mesh = slim_cylinder(length, diameter_base / 2.0, diameter_top / 2.0)

    return mesh


def compute_element(element_node, leaves, min_length=0.01, classic=False):
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None

    if n.label.startswith("Leaf"):  # leaf element
        blade = n.complex()
        species = blade.species
        if blade.visible_length >= min_length:
            if blade.shape_key is not None and n.srb is not None:
                if leaves[species].dynamic:
                    inclin = 1  # inclination is encoded in db
                else:
                    inclin = blade.inclination
                geom = leaves[
                    species
                ].mesh(
                    blade.shape_key,
                    blade.shape_mature_length,
                    blade.shape_max_width,
                    blade.visible_length,
                    n.srb,
                    n.srt,
                    incline=inclin,
                    flipx=True,
                    min_area=min_length**2 / 100,
                )  # flipx allows x-> -x to place the shape along with the tiller positioned with turtle.down()

            if n.lrolled > 0:
                rolled = StemElement_mesh(n.lrolled, n.d_rolled, n.d_rolled, classic)
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        stem = n.complex()
        # diameter_base = stem.parent().diameter if (stem.parent() and stem.parent().diameter > 0.) else stem.diameter
        # diameter_top = n.diam
        if n.length >= min_length:
            diameter_base = stem.diameter
            diameter_top = stem.diameter
            geom = StemElement_mesh(n.length, diameter_base, diameter_top, classic)

    return geom


class AdelTurtle(pgl.PglTurtle):
    def __init__(self):
        super(AdelTurtle, self).__init__()
        self.context = {}

    def transform(self, mesh, face_up=False):
        x = self.getUp()
        if face_up:
            z = pgl.Vector3(0, 0, 1)
        else:
            z = self.getHeading()
        bo = pgl.BaseOrientation(x, z ^ x)
        matrix = pgl.Transform4(bo.getMatrix())
        matrix.translate(self.getPosition())
        # print 'Position ', turtle.getPosition()
        mesh = mesh.transform(matrix)
        return mesh

    def getFrame(self):
        pos = self.getPosition()
        Head = self.getHeading()
        Up = self.getUp()
        return {"Position": pos, "Head": Head, "Up": Up}

    def setFrame(self, frame):
        self.move(frame["Position"])
        self.setHead(frame["Head"], frame["Up"])


class AdelVisitor:
    """Performs geometric interpretation of mtg nodes"""

    def __init__(self, leaves, min_length, classic, face_up):
        self.classic = classic
        self.face_up = face_up
        self.min_length = min_length
        self.leaves = leaves

    def __call__(self, g, v, turtle):
        geometry = g.property("geometry")
        # 1. retrieve the node
        n = g.node(v)
        axis = n.complex_at_scale(scale=2).label
        az_axis = n.complex_at_scale(scale=2).azimuth
        prev_axis = turtle.context.get("axis", axis)
        metamer = n.complex_at_scale(scale=3)

        # Go to plant position if first plant element
        if n.parent() is None:  # this is a new plant base
            p = n.complex_at_scale(scale=1)
            if "position" in p.properties():
                # print p.label, 'moving to ', p.position
                turtle.move(list(map(float, p.position)))
            else:
                turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)
            if "azimuth" in p.properties():
                turtle.rollR(p.azimuth)
            # print prev_axis, 'stop'
            # print 'new MS start'
            turtle.context.update(
                {
                    "MS_top": turtle.getFrame(),
                    "tiller_base": turtle.getFrame(),
                    "top": turtle.getFrame(),
                    "is_axis_first_StemElement": True,
                }
            )

        if axis != prev_axis:
            if metamer.edge_type() == "+":
                turtle.context["is_axis_first_StemElement"] = True
                if prev_axis == "MS":
                    # this is the begining of the first tiller
                    # print axis, 'start'
                    # top of mainstem saved
                    turtle.context["MS_top"] = turtle.context["top"]
                    turtle.context["tiller_base"] = turtle.context["top"]
                else:
                    # this is a new tiller attached to the same point thatn the previous tiller
                    # print prev_axis, 'stop'
                    # print axis, 'start'
                    # return to tilller base
                    turtle.context["top"] = turtle.context["tiller_base"]
            else:  # this is the continuation of MS
                # print prev_axis, 'stop'
                # print axis, 'continue'
                turtle.context["top"] = turtle.context["MS_top"]

        # hypothesis that inclin is to be applied at the base of the visible elements
        # if n.offset > 0:
        #    turtle.f(n.offset)
        turtle.setFrame(turtle.context["top"])
        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            inclin = float(n.inclination) if n.inclination else 0.0
            azim = float(n.azimuth) if n.azimuth else 0.0
            if inclin:
                # print 'node', n._vid, 'inclin', inclin
                # incline along curent azimuth for ramification (tiller bases) or plant base
                if turtle.context["is_axis_first_StemElement"]:
                    # print 'axis',axis, 'prev_axis', prev_axis,' node ', n._vid, 'edge', n.edge_type(),'up before inclin ', turtle.getUp(), 'inclin', inclin
                    if (
                        prev_axis != "MS"
                    ):  # new tiller attached to the same position than the firt
                        turtle.rollR(az_axis)
                        turtle.context["tiller_base"] = turtle.getFrame()
                    turtle.down(inclin)
                    turtle.context["is_axis_first_StemElement"] = False
                    # print 'up after inclin', turtle.getUp()
                # if not incline towardss vertical
                else:
                    up = turtle.getUp()
                    zleft = turtle.getLeft()[2]
                    turtle.rollToVert()
                    # print 'up after rollToVert', turtle.getUp()
                    angle = degrees(pgl.angle(up, turtle.getUp()))
                    dzl = zleft - turtle.getLeft()[2]
                    turtle.down(inclin)
                    # replace turtle in original azimuth plane
                    # print 'angle ', angle, 'dzl', dzl
                    if dzl < 0:
                        angle = -angle
                    turtle.rollR(-angle)
            if azim:
                # print 'node', n._vid, 'azim ', azim
                turtle.rollR(azim)

        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            # update geometry of elements
            mesh = None
            if n.length > 0:
                mesh = compute_element(n, self.leaves, min_length=self.min_length, classic=self.classic)
            if mesh:
                n.geometry = turtle.transform(
                    mesh, face_up=self.face_up and n.label.startswith("Leaf")
                )
                n.anchor_point = turtle.getPosition()
            else:
                if v in geometry:  # delete existing geometry
                    geometry.pop(v)
        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.length > 0:
                turtle.f(n.length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})
        turtle.context.update({"axis": axis})


def mtg_interpreter(g, leaves, min_length=0.01, classic=False, face_up=False):
    """Compute/update the geometry on each node of the MTG using Turtle geometry."""
    # BUG : sub_mtg mange le vertex plant => on perd la plante !
    # plants = g.component_roots_at_scale(g.root, scale=1)
    # nplants = g.nb_vertices(scale=1)
    # gt = MTG()

    # for plant in plants:
    #   gplant = g.sub_mtg(plant)
    turtle = AdelTurtle()
    visitor = AdelVisitor(leaves, min_length, classic, face_up)
    scene = TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False, all_roots=True)
    #   gt = union(gplant,gt)

    return g


def plot3d(g, leaf_material=None, stem_material=None, soil_material=None, colors=None):
    """
    Returns a plantgl scene from an mtg.
    """

    Material = pgl.Material
    Color3 = pgl.Color3
    Shape = pgl.Shape
    Scene = pgl.Scene

    if colors is None:
        if leaf_material is None:
            leaf_material = Material(Color3(0, 180, 0))
        if stem_material is None:
            stem_material = Material(Color3(0, 130, 0))
        if soil_material is None:
            soil_material = Material(Color3(170, 85, 0))
        # colors = g.property('color')

    geometries = g.property("geometry")
    greeness = g.property("is_green")
    labels = g.property("label")
    scene = Scene()

    def geom2shape(vid, mesh, scene, colors):
        shape = None
        if isinstance(mesh, list):
            for m in mesh:
                geom2shape(vid, m, scene, colors)
            return
        if mesh is None:
            return
        if isinstance(mesh, Shape):
            shape = mesh
            mesh = mesh.geometry
        label = labels.get(vid)
        is_green = greeness.get(vid)
        if colors:
            shape = Shape(mesh, Material(Color3(*colors.get(vid, [0, 0, 0]))))
        elif not greeness:
            if not shape:
                shape = Shape(mesh)
        elif label.startswith("Stem") and is_green:
            shape = Shape(mesh, stem_material)
        elif label.startswith("Leaf") and is_green:
            shape = Shape(mesh, leaf_material)
        elif not is_green:
            shape = Shape(mesh, soil_material)
        shape.id = vid
        scene.add(shape)

    for vid, mesh in geometries.items():
        geom2shape(vid, mesh, scene, colors)
    return scene


def transform_geom(geom, translation, rotation):
    # force cast to float (pgl does not accept values extracted from numpy arryas
    translation = list(map(float, translation))
    if isinstance(geom, pgl.Geometry):
        geom = pgl.Translated(translation, pgl.AxisRotated((0, 0, 1), rotation, geom))
    elif isinstance(geom, pgl.Shape):
        geom = pgl.Shape(
            pgl.Translated(
                translation, pgl.AxisRotated((0, 0, 1), rotation, geom.geometry)
            )
        )
    return geom
