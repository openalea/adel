from rpy import r
import os
from os.path import join
from openalea.core.pkgmanager import PackageManager

pm = PackageManager()
pkg = pm.get("openalea.adel.adelMons")

path_ini = os.getcwd()

path = ""
if pkg:
    path = pkg.path


class write_graphtalPars:
    """write a .h input file for adelSimexp from excel files"""

    def __init__(self):
        pass

    def __call__(
        self, gen, semis, morpho_f, phyllo_f, sen_f, prev_f, alpha_f, plan_f, fname
    ):
        os.chdir(path)
        r.source("writeparsGL.R")
        r.readExcell.local_mode(0)
        r.as_data_frame.local_mode(0)

        bloc = "B1"  # B1 par defaut (changer si bcp de variabilite entre blocs)

        # dimensions des organes (en, gaine, longueur largeur de feuilles) f(gen,semis)
        morpho = r.readExcell(morpho_f, gen + "_" + semis)

        # tableau angles des feuilles  f(gen,semis)
        anglesPrev = r.readExcell(prev_f, gen + "_" + semis)

        # dynamique nbvis et nblig
        dev = r.readExcell(phyllo_f, gen + "_" + semis)
        n_max = round(dev.as_py()["nbvis"][-1])

        # dynamique de senescence
        r.readExcell.local_mode(1)
        sen = r.readExcell(sen_f, gen + "_" + semis)
        r.readExcell.local_mode(0)
        sen["TT"].append(3000.0)
        sen["nb_sen"].append(n_max + 1)
        sen = r.as_data_frame(sen)

        # formes de limbes
        alpha = r.readExcell(alpha_f, gen)

        # plan parcelle
        plan = r.readExcell(plan_f, gen + "_" + semis + "_" + bloc)

        # fname = "_05_geomF2_S1B1.h"

        txt1 = r.writenumf(fname, dev, morpho, anglesPrev)
        print(txt1)
        r.writePO(fname)
        r.writecarto(fname, plan, 0.8, 0.125)
        txt2 = r.writedimkin(fname, dev, sen, morpho)
        print(txt2)
        r.writeprev(fname, alpha, anglesPrev, dev)

        os.chdir(path_ini)

        return join(join(path, "temp"), fname)
