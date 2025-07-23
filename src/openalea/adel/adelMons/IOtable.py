def N_lignes(fichier):
    """compte le nombre de lignes d'un fichier (compte le nombre d'elements de la liste readlines()"""

    position_ini = fichier.tell()
    fichier.seek(0)
    N_tot_lignes = len(fichier.readlines())
    fichier.seek(position_ini)
    return N_tot_lignes


def transcript_csv_str(fichier):
    """transcrit une ligne d'un fichier.csv en liste de float"""
    ligne_ch = fichier.readline()
    liste, i, caract = (
        [],
        0,
        "",
    )  # initialise une liste vide et un indice a 0 et une chaine de caractere vide

    if ligne_ch == "":  # verifie si la chaine est vide
        return "chaine vide"

    while ligne_ch[i] != "\n":  # boucle pour une valeur en .csv (separateur ;)
        ch = ligne_ch[i]
        if ligne_ch[i] == ",":  # converti separateur decimal en point
            ch = "."
        if ligne_ch[i] == ";":
            liste.append(caract)
            caract = ""
        if ligne_ch[i] != ";":
            caract = caract + ch
        i = i + 1

    liste.append(caract)  # ajoute la derniere valeur
    return liste


def table_csv_str(fichier):
    """transcrit un fichier .csv de nombres en un tableau (liste de liste) de nombre reels"""

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)

    for i in range(n):
        liste.append(transcript_csv_str(fichier))

    return liste


def transcript_csv(fichier):
    """transcrit une ligne d'un fichier.csv en liste de float"""
    ligne_ch = fichier.readline()
    liste, i, caract = (
        [],
        0,
        "",
    )  # initialise une liste vide et un indice a 0 et une chaine de caractere vide

    if ligne_ch == "":  # verifie si la chaine est vide
        return "chaine vide"

    while ligne_ch[i] != "\n":  # boucle pour une valeur en .csv (separateur ;)
        ch = ligne_ch[i]
        if ligne_ch[i] == ",":  # converti separateur decimal en point
            ch = "."
        if ligne_ch[i] == ";":
            liste.append(float(caract))
            caract = ""
        if ligne_ch[i] != ";":
            caract = caract + ch
        i = i + 1

    liste.append(float(caract))  # ajoute la derniere valeur

    return liste


def table_csv(fichier):
    """transcrit un fichier .csv de nombres en un tableau (liste de liste) de nombre reels"""

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)

    for i in range(n):
        liste.append(transcript_csv(fichier))

    return liste


def transcript_txt(fichier):
    """transcrit une ligne d'un fichier.txt en liste de chaine"""
    ligne_ch = fichier.readline()
    liste, i, caract = (
        [],
        0,
        "",
    )  # initialise une liste vide et un indice a 0 et une chaine de caractere vide

    if ligne_ch == "":  # verifie si la chaine est vide
        return "chaine vide"

    while ligne_ch[i] != "\n":  # boucle pour une valeur en .csv (separateur ;)
        ch = ligne_ch[i]
        if (ligne_ch[i] == " " or ligne_ch[i] == "\t") and caract != "":
            liste.append(caract)
            caract = ""
        if ligne_ch[i] != " " and ligne_ch[i] != "\t":
            caract = caract + ch
        i = i + 1

    liste.append(caract)  # ajoute la derniere valeur
    return liste


def table_txt(fichier):
    """transcrit un fichier .txt un tableau (liste de liste) de chaines"""

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)

    for i in range(n):
        liste.append(transcript_txt(fichier))

    return liste


def ecriture_csv(table, fichier):
    """ecrit table (liste de liste) de chiffre dans un fichier csv"""
    for i in range(len(table)):
        for j in range(len(table[i][0:]) - 1):
            fichier.write(str(table[i][j]))
            fichier.write(";")
        fichier.write(str(table[i][j + 1]))
        fichier.write("\n")
    fichier.close()


def ecriture_txt(table, fichier):
    """ecrit table (liste de liste) de chiffre dans un fichier txt"""
    for i in range(len(table)):
        for j in range(len(table[i][0:]) - 1):
            fichier.write(str(table[i][j]))
            fichier.write(" ")
        fichier.write(str(table[i][j + 1]))
        fichier.write("\n")
    fichier.close()


def copie_partielle(fichier, out, n, m):
    """recopie de fichier des lignes n a m dans le fichier out"""
    for i in range(m):
        ligne_ch = fichier.readline()
        if i >= n:
            out.writelines(ligne_ch)
