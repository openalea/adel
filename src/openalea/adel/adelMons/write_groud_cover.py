from os.path import join

from openalea.core.pkgmanager import PackageManager

from . import IOtable

pm = PackageManager()
pkg = pm.get("openalea.adel.adelMons")
path = pkg.path


def write_groud_cover(table, file_name):
    """write ground cover results"""
    res = []
    for i in range(len(table)):
        for j in range(len(table[i])):
            res.append(table[i][j])

    # ecriture du fichier
    name = join(path, file_name)
    out = file(name, "w")
    IOtable.ecriture_csv(res, out)
    out.close()

    return name
