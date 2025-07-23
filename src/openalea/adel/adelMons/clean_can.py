from . import IOtable
from scipy import sqrt


def distance(p1, p2):
    return sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def triangle_area(p1, p2, p3):
    """formule aire triangle a partir des coord de ses sommets"""
    a = distance(p1, p2)
    b = distance(p1, p3)
    c = distance(p2, p3)
    return 0.25 * sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c))


class clean_can:
    """remove abnormal triangles genetared in graphtal by the 'atan' bug"""

    def __init__(self):
        pass

    def __call__(self, name, surf_lim, hlim):
        f = file(name, "r")
        tab_geom = IOtable.table_txt(f)
        f.close()
        for i in range(1, len(tab_geom) - 1):
            tab_geom[i][5:] = list(map(float, tab_geom[i][5:]))

        out = []
        out.append(tab_geom[0])
        for i in range(1, len(tab_geom) - 1):
            if (
                triangle_area(tab_geom[i][5:8], tab_geom[i][8:11], tab_geom[i][11:])
                < surf_lim
                and abs(tab_geom[i][7]) < hlim
                and abs(tab_geom[i][10]) < hlim
                and abs(tab_geom[i][13]) < hlim
            ):  # test pour ecarter triangles bugge
                out.append(tab_geom[i])

        out.append(tab_geom[-1])

        f = file(name, "w")
        IOtable.ecriture_txt(out, f)
        f.close()

        return name
