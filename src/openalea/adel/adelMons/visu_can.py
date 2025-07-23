from . import IOtable
from openalea.plantgl.all import *


class visu_can:
    """can file visualisation"""

    def __init__(self):
        pass

    def __call__(self, can_file):
        # recup du fichier can dans tb
        f = file(can_file, "r")
        tb = IOtable.table_txt(f)
        f.close()

        # retire premiere et dernieres lignes
        tb = tb[1:-1]

        # cree une scene vide, un viewer
        MaScene = Scene()
        MonViewer = Viewer

        # ajout des triangles 2 par 2 pour eviter bugs lors de l'exportation des triangleset dans pov-ray
        for i in range(0, len(tb) - 1, 2):
            # creation d'une liste de coordonnees des points et d'une liste d'index des triangles pour les lignes i et i+1
            ind, pts = [], []
            # count=0

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
            if tb[i][2][0] != "2":  # sp_opt feuille senescente
                R, G, B = 0, 160, 0  # vert
            else:
                R, G, B = 255, 204, 0  # jaune

            MaScene.add(Shape(TriangleSet(points, indices), Material(Color3(R, G, B))))

        # Affichage de la scene
        MonViewer.display(MaScene)

        return can_file
