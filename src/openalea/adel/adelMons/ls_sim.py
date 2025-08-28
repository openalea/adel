import math
import os

from openalea.core.pkgmanager import PackageManager

from . import IOtable

pm = PackageManager()
pkg = pm.get("openalea.adel.adelMons")
path = ""
if pkg:
    path = pkg.path

path_ini = os.getcwd()


def t_list(tab):
    """transpose tab"""
    res = []
    for j in range(len(tab[0])):
        v = []
        for i in range(len(tab)):
            v.append(tab[i][j])

        res.append(v)

    return res


def conv_dataframe(tab):
    """converti liste de liste en dictionnaire; prend la cle comme le pemier element de la liste"""
    dat = {}
    for i in range(len(tab)):
        dat[str(tab[i][0])] = tab[i][1:]

    return dat


def extract_dataframe(dat, ls_cles, cle, val=None):
    """extrait dans listes de cles ls_cles les lignes pour lesquelles cle=val; toutes si val=None"""
    # cree liste d'index ou cle = val

    id = []
    for i in range(len(dat[cle])):
        if val is None:
            id.append(i)
        else:
            if dat[cle][i] == val:
                id.append(i)

    x = {}
    for k in ls_cles:  # recupere les paires interessantes
        v = []
        for i in id:  # les id respectant cle=val
            v.append(dat[k][i])

        x[k] = v

    return x


class ls_sim:
    """Gets the lists of argument to pass for batch runs with adelsimexp.lsys (list of TT, list of can file names)
    Converts fov angle into camera angle for povray (enables to simulate ground cover with povray)"""

    def __init__(self):
        pass

    def __call__(self, ls_file):
        f = open(ls_file, "r")
        ls = IOtable.table_csv_str(f)
        f.close()
        ls = conv_dataframe(t_list(ls))
        ls = extract_dataframe(ls, list(ls.keys()), "simul", val="1")
        ls["TT"] = list(map(int, ls["TT"]))

        # generer noms de fichier param / noms de fichiers can
        ls_par, ls_can, ls_pos, ls_angle = [], [], [], []
        # pourrait etre passe en parametre de la boite ls_sim; de meme pour le calcul ou non de l'angle d'ouverture pov-ray

        os.chdir(path)
        for i in range(len(ls["TT"])):
            ls_par.append(
                ".\GenPars"
                + "\\"
                + ls["annee"][i][-2:]
                + "_geom"
                + ls["gen"][i]
                + "_"
                + ls["semis"][i]
                + ls["bloc"][i]
                + ".h"
            )
            ls_can.append(
                "temp"
                + "\\"
                + ls["gen"][i]
                + "_"
                + ls["semis"][i]
                + "_"
                + ls["annee"][i][-2:]
                + "_TT"
                + str(ls["TT"][i])
                + ".can"
            )
            ls_pos.append(float(ls["Z camera"][i]))
            angle = 2.0 * math.degrees(
                math.atan2(
                    0.8 * math.tan(math.radians(float(ls["angle ouverture"][i])) / 2.0),
                    1.0,
                )
            )  # angle d'ouverture de la camera pov-ray correspondante
            ls_angle.append(angle)

        os.chdir(path_ini)

        return ls_par, ls["TT"], ls_can, ls_angle, ls_pos
