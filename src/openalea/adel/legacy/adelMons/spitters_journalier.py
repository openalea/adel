from math import sin, cos, radians


def DecliSun(JJ):
    """Declinaison (rad) du soleil en fonction du jour julien"""
    alpha = 2 * 3.14 * (JJ - 1) / 365
    return 0.006918 - 0.399912 * cos(alpha) + 0.070257 * sin(alpha)


def extra(Rg, JJ, heureTU, latitude):
    """rayonnement extraterrestre horarire"""
    hrad = 2 * 3.14 / 24 * (heureTU - 12)
    lat = radians(latitude)
    dec = DecliSun(JJ)
    costheta = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrad)
    Io = (
        1370 * (1 + 0.033 * cos(2 * 3.14 * (JJ - 4) / 366))
    )  # eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
    So = Io * costheta  # eclairement dans un plan parallele a la surface du sol
    return So


def spitters_journalier(Rg, JJ, latitude):
    """fraction diffus/Global en fonction du rapport Global(Sgd)/Extraterrestre(Sod)- pas de temps journalier"""
    """ Rg et So en J.m-2.d-1 """
    lat = radians(latitude)
    dec = DecliSun(JJ)
    So = 0
    for heureTU in range(24):
        hrad = 2 * 3.14 / 24 * (heureTU - 12)
        costheta = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrad)
        Io = (
            1370 * (1 + 0.033 * cos(2 * 3.14 * (JJ - 4) / 366))
        )  # eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
        So = So + max(
            0.0, Io * costheta * 3600
        )  # eclairement dans un plan parallele a la surface du sol, 0. la nuit

    RsRso = Rg / So
    if RsRso < 0.07:
        return 1
    elif (RsRso >= 0.07) and (RsRso < 0.35):
        return 1 - 2.3 * (RsRso - 0.07) * (RsRso - 0.07)
    elif (RsRso >= 0.35) and (RsRso < 0.75):
        return 1.33 - 1.46 * RsRso
    else:
        return 0.23
