from . import IOtable
from os.path import join


def t_list(tab):
    """transpose tab"""
    res = []
    for j in range(len(tab[0])):
        v = []
        for i in range(len(tab)):
            v.append(tab[i][j])

        res.append(v)

    return res


def conv_list(tab):
    """converti dictionnaireen liste de liste en ;  cle comme pemier element de la liste"""
    """ format compatible pour mes_csv"""
    dat = []
    for i in list(tab.keys()):
        v = [i]
        dat.append(v)

    count = 0
    for i in list(tab.keys()):
        for j in range(len(tab[i])):
            dat[count].append(tab[i][j])

        count = count + 1

    return dat


def conv_list2(tab):
    """converti dictionnaire avec 1 seule valeur par cle en liste de liste ;  cle comme pemier element de la liste"""
    """ format compatible pour mes_csv"""
    dat = []
    for i in list(tab.keys()):
        v = [i, tab[i]]
        dat.append(v)

    return dat


def write_dict(dict, directory, name):
    try:
        tab = t_list(conv_list(dict))
    except:
        tab = conv_list2(dict)

    out = file(join(directory, name), "w")
    IOtable.ecriture_csv(tab, out)
    out.close()

    return join(directory, name)
