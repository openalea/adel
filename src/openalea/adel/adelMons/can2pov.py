from . import IOtable
from openalea.plantgl.all import *


class can2pov:
    """converts can file into pov file using geom viewer"""

    def __init__(self):
        pass

    def __call__(self, can_file, cam_type, backg, soil, cam_pos, fov, cam_rot):
        # recup du fichier can dans tb
        f = file(can_file, "r")
        tb = IOtable.table_txt(f)
        f.close()

        # retire premiere et dernieres lignes
        tb = tb[1:-1]

        # cree une scene vide, un viewer
        MaScene = Scene()

        # ajout des triangles 2 par 2 pour eviter bugs lors de l'exportation des triangleset dans pov-ray
        for i in range(0, len(tb) - 1, 2):
            indices = Index3Array([(0, 1, 2), (3, 4, 5)])
            coord1 = list(map(float, tb[i][5:]))
            coord2 = list(map(float, tb[i + 1][5:]))
            points = Point3Array(
                [
                    Vector3(*coord1[:3]),
                    Vector3(*coord1[3:6]),
                    Vector3(*coord1[6:9]),
                    Vector3(*coord2[:3]),
                    Vector3(*coord2[3:6]),
                    Vector3(*coord2[6:9]),
                ]
            )
            if tb[i][3] != "8":  # 8 pour limbes senescents et panicule
                R, G, B = 0, 160, 0  # vert
            else:
                R, G, B = 255, 204, 0  # jaune

            MaScene.add(Shape(TriangleSet(points, indices), Material(Color3(R, G, B))))

        # exporte scene en format pov:
        # cree fichier qui contient les mesh
        t = Tesselator()
        pov = PovFilePrinter(can_file[0:-4] + "_mesh.pov", t)
        MaScene.apply(pov)

        # cree fichier pov devinitif
        f = file(can_file[0:-4] + ".pov", "w")
        self.pov_header(f, cam_type, backg, soil, cam_pos, fov, cam_rot)
        f.write('#include "' + can_file[0:-4] + "_mesh.pov" + '"\n')
        f.close()

        return can_file[0:-4] + ".pov"

    @staticmethod
    def pov_header(f, cam_type, backg, soil, cam_pos, fov, cam_rot):
        """cree une camera povray, un arriere plan et lumiere"""
        """ sol optionnel """

        f.write(
            "/////////////////////////////////////////////////////////////////////////////////"
            + "\n"
        )
        f.write('#include "colors.inc"' + "\n\n")

        f.write("camera {" + "\n")
        f.write("   " + cam_type + "\n")
        f.write(
            "    location <"
            + str(cam_pos[0])
            + ","
            + str(cam_pos[1])
            + ","
            + str(cam_pos[2])
            + ">"
            + "\n"
        )
        f.write("    direction <-1,0,0>" + "\n")
        f.write("    up <0,0,1>" + "\n")
        f.write("    right <0,4/3,0>" + "\n")
        f.write("    angle " + str(fov) + "    " + "\n")
        f.write(
            "    rotate <"
            + str(cam_rot[0])
            + ","
            + str(cam_rot[1])
            + ","
            + str(cam_rot[2])
            + ">}"
            + "\n\n"
        )  # <0,-16,-89> pour film, par defaut; <0,-90,90> pour taux de couv

        f.write("light_source {  <10.6066,0,10.6066> color rgb <1,1,1>}" + "\n\n")

        f.write(
            "background { color rgb <"
            + str(backg[0])
            + ","
            + str(backg[1])
            + ","
            + str(backg[2])
            + "> }"
            + "\n\n"
        )

        if soil:
            f.write("#declare T1=" + "\n")
            f.write("    texture {" + "\n")
            f.write("      pigment { color red 1.0 green 0.75 blue 0.33 }" + "\n")
            f.write("      normal {bumps 0.4  scale 0.007}}" + "\n\n")

            f.write("#declare T2=" + "\n")
            f.write(" texture{" + "\n")
            f.write("   pigment{color red 1.0 green 0.75 blue 0.33}" + "\n")
            f.write("   normal{granite 0.6 scale 0.1 }" + "\n")
            f.write("   finish{phong 0.8 phong_size 200}}   " + "\n\n")

            f.write("box {" + "\n")
            f.write("    <-3, -2,   0>,  " + "\n")
            f.write("    < 0.5, 0.5,  0.001>   " + "\n")
            f.write("    texture {T2}}" + "\n\n")

            # inclure light source parameters?
