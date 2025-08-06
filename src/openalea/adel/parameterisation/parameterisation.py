from openalea.adel import AdelR
from openalea.adel.fitting import fit2

import numpy
from rpy2 import robjects

from rpy2.robjects import numpy2ri

# try:
#     numpy2ri.activate()  # force auto-conversion mode of Robject to array
# except:
#     pass


class PlantParameter:
    def __init__(self, **kwargs):
        pass


def dataframe(d):
    df = {}
    if d is None:
        return robjects.r("as.null()")
    else:
        for k, v in d.items():
            df[k] = numpy2ri(numpy.array(v))
    dataf = robjects.r["data.frame"](**df)
    return dataf


def simpleMais(simpleMais_parameter):
    """generate Adel canopy Table (R)  from a simpleMais parameter dictionary"""

    d = simpleMais_parameter
    nrow = len(d["plante"])
    zero = numpy.zeros(nrow)

    if ("longEn" not in d) or ("longGa" not in d):
        d.update({"longEn": zero})
        d.update({"longGa": (numpy.array(d["longTige"])).cumsum()})

    canT = {
        "plant": d["plante"],
        "axe": zero,
        "numphy": d["phytomere"],
        "Lv": d["longFeu"],
        "Ll": d["longFeu"],
        "L_shape": d["longFeu"],
        "Lsen": zero,
        "Lw": d["largFeu"],
        "Lw_shape": d["largFeu"],
        "LcType": d["phytomere"],
        "LcIndex": zero,
        "Linc": d["inclinaisonFeu"],
        "Laz": d["azimuthFeu"],
        "Lpo": zero + 1,
        "Lpos": zero + 2,
        "Gl": d["longGa"],
        "Gv": d["longTige"],
        "Gsen": zero,
        "Gpo": zero + 1,
        "Gpos": zero + 2,
        "Gd": d["diamTige"],
        "Ginc": zero,
        "El": d["longEn"],
        "Ev": zero,
        "Esen": zero,
        "Ed": zero,
        "Einc": zero,
        "Epo": zero + 1,
        "Epos": zero + 2,
    }
    return canT
    # return(dataframe(canT))


def simpleMais2dict(simpleMais_parameter):
    """generate Adel canopy Table (R)  from a simpleMais parameter dictionary"""

    d = simpleMais_parameter
    nrow = len(d["plante"])
    zero = numpy.zeros(nrow)

    if ("longEn" not in d) or ("longGa" not in d):
        d.update({"longEn": zero})
        d.update({"longGa": (numpy.array(d["longTige"])).cumsum()})

    canT = {
        "plant": d["plante"],
        "axe": zero,
        "numphy": d["phytomere"],
        "Lv": d["longFeu"],
        "Ll": d["longFeu"],
        "L_shape": d["longFeu"],
        "Lsen": zero,
        "Lw": d["largFeu"],
        "Lw_shape": d["largFeu"],
        "LcType": d["phytomere"],
        "LcIndex": zero,
        "Linc": d["inclinaisonFeu"],
        "Laz": d["azimuthFeu"],
        "Lpo": zero + 1,
        "Lpos": zero + 2,
        "Gl": d["longGa"],
        "Gv": d["longTige"],
        "Gsen": zero,
        "Gpo": zero + 1,
        "Gpos": zero + 2,
        "Gd": d["diamTige"],
        "Ginc": zero,
        "El": d["longEn"],
        "Ev": zero,
        "Esen": zero,
        "Ed": zero,
        "Einc": zero,
        "Epo": zero + 1,
        "Epos": zero + 2,
    }
    return (canT,)


def bell_shaped_dist(total_area=1, nb_phy=15, rmax=0.7, skewness=5):
    """returns leaf area of individual leaves along bell shaped model"""

    r = numpy.linspace(1.0 / nb_phy, 1, nb_phy)
    k = skewness
    relative_surface = numpy.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    leaf_area = relative_surface / relative_surface.sum() * total_area
    return leaf_area.tolist()


def shape_factor(leaf_database, nb_phy):
    """compute shape factors for nb_phy in a leaf database"""

    db = leaf_database
    rank_max = max(int(x) for x in db)

    def choose_rank(rank):
        if str(rank) not in db:
            rank = rank + 1 if str(rank + 1) in db else rank - 1
        if str(rank) in db:
            return rank

        while str(rank) not in db and rank < rank_max:
            rank += 1
        if str(rank) not in db:
            rank = rank_max
        return rank

    ranks = [choose_rank(n) for n in range(1, nb_phy + 1)]
    norm_surface = numpy.array([fit2(*db[str(rank)][0])[1] for rank in ranks])

    return norm_surface


def geometric_dist(height=15, nb_phy=15, q=1):
    """returns distances between individual leaves along a geometric model"""

    if q == 1:
        u0 = float(height) / nb_phy
    else:
        u0 = height * (1.0 - q) / (1.0 - q ** (nb_phy + 1))

    return [u0 * q**i for i in range(nb_phy)]


def simpleMais_param(
    total_area=10000,
    total_height=200,
    pseudostem_height=20,
    nb_phy=16,
    nb_young_phy=6,
    # leaf area dist function
    lad_skew=5,
    lad_rmax=0.67,
    # proression of distances between leaves
    pseudostem_dist=1.4,
    stem_dist=1.0,
    # Phyllotaxy for young and mature phytomer
    phyllotactic_angle=180,
    # insertion angle
    basal_insertion=50,
    top_insertion=30,
    diam_base=2.5,
    diam_top=1,
    leaf_width_ratio=0.1,
    shape_factors=0.75,
):
    """generate simpleMais dictionary from parameter"""

    nb_phy = int(nb_phy)
    nb_young_phy = int(nb_young_phy)

    # compute the leaf surface
    leaf_area = numpy.array(bell_shaped_dist(total_area, nb_phy, lad_rmax, lad_skew))
    # derive corresponding length and widths
    lengths = numpy.sqrt(leaf_area / numpy.array(shape_factors) / leaf_width_ratio)
    widths = lengths * leaf_width_ratio

    # distances between leaves
    pseudostem = geometric_dist(pseudostem_height, nb_young_phy, pseudostem_dist)
    stem = geometric_dist(
        total_height - pseudostem_height, nb_phy - nb_young_phy, stem_dist
    )
    distances = pseudostem + stem
    # sheath and internode true length
    pseudostem_sheath = (numpy.array(pseudostem)).cumsum()
    pseudostem_internode = [0] * nb_young_phy
    stem_sheath = numpy.array(stem)
    stem_internode = [0] + (stem_sheath[1:]).tolist()
    sheaths = pseudostem_sheath.tolist() + stem_sheath.tolist()
    internodes = pseudostem_internode + stem_internode
    # stem diameters
    diameters = [diam_base] * nb_young_phy + numpy.linspace(
        diam_base, diam_top, nb_phy - nb_young_phy
    ).tolist()
    # inclinations
    inclinations = numpy.linspace(basal_insertion, top_insertion, nb_phy)
    # azimuths
    az0 = 0
    azimuths = [az0 + i * phyllotactic_angle for i in range(nb_phy)]
    # output dict
    dTags = [
        "plante",
        "phytomere",
        "longFeu",
        "largFeu",
        "inclinaisonFeu",
        "azimuthFeu",
        "longTige",
        "diamTige",
        "longGa",
        "longEn",
    ]
    dVals = [
        [1] * nb_phy,
        list(range(1, nb_phy + 1)),
        lengths,
        widths,
        inclinations,
        azimuths,
        distances,
        diameters,
        sheaths,
        internodes,
    ]
    dVals = list(map(numpy.array, dVals))
    dout = dict(list(zip(dTags, dVals)))

    return dout


def MonoAxeWheat_param(
        axedim=None,
    inclination=30,
    scale_stem=1,
    scale_leaf=1,
    scale_stem_diameter=1,
    scale_leaf_width=1,
):
    """
    Generate a simple wheat axe geometric model with leaves and stems scaled from an Axedim dictionary and fixed azimutal parameterisation
    """

    if axedim is None:
        axedim = {
            "Lamina_length": [
                8.125,
                9.25,
                9.35,
                10,
                11.4,
                13.7,
                16.55,
                19.8,
                25.175,
                28.8,
                24.1,
            ],
            "Lamina_width": [0.3, 0.325, 0.4, 0.45, 0.55, 0.75, 1, 1.2, 1.28,
                             1.425, 1.8],
            "Sheath_length": [
                3,
                3.05,
                3.05,
                3.4,
                4.2,
                6.225,
                9.125,
                12,
                14.2,
                17.2,
                18.675,
            ],
            "Internode_length": [0, 0, 0, 0, 0, 0.1, 1.9, 6.1, 9.675, 14.45,
                                 16.95],
            "Stem_diameter": [
                0.14,
                0.18,
                0.21,
                0.24,
                0.29,
                0.34,
                0.36,
                0.4,
                0.48,
                0.54,
                0.73,
            ],
        }
    nb_phy = len(axedim["Lamina_length"])

    hfeu = numpy.array([0] + axedim["Sheath_length"]) + numpy.array(
        [0] + axedim["Internode_length"]
    )

    dTags = [
        "plante",
        "phytomere",
        "longFeu",
        "largFeu",
        "inclinaisonFeu",
        "azimuthFeu",
        "longTige",
        "diamTige",
    ]
    dVals = [
        [1] * nb_phy,
        list(range(1, nb_phy + 1)),
        numpy.array(axedim["Lamina_length"]) * scale_leaf,
        numpy.array(axedim["Lamina_width"]) * scale_leaf_width,
        [inclination] * nb_phy,
        [i * 180 for i in range(nb_phy)],
        numpy.diff(hfeu) * scale_stem,
        numpy.array(axedim["Stem_diameter"]) * scale_stem_diameter,
    ]

    dVals = list(map(numpy.array, dVals))
    d = dict(list(zip(dTags, dVals)))

    nrow = len(d["plante"])
    zero = numpy.zeros(nrow)

    canT = {
        "plant": d["plante"],
        "axe": zero,
        "numphy": d["phytomere"],
        "Lv": d["longFeu"],
        "Ll": d["longFeu"],
        "Lsen": zero,
        "Lw": d["largFeu"],
        "LcType": zero + 1,
        "LcIndex": zero,
        "Linc": d["inclinaisonFeu"],
        "Laz": d["azimuthFeu"],
        "Lpo": zero + 1,
        "Lpos": zero + 2,
        "Gl": zero,
        "Gv": zero,
        "Gsen": zero,
        "Gpo": zero + 1,
        "Gpos": zero + 2,
        "Gd": zero,
        "Ginc": zero,
        "El": d["longTige"],
        "Ev": d["longTige"],
        "Esen": zero,
        "Ed": d["diamTige"],
        "Einc": zero,
        "Epo": zero + 1,
        "Epos": zero + 2,
    }
    return canT


def simpleMais_plan(d, step=0):
    params = [
        "S",
        "H",
        "Hps",
        "NbPhy",
        "NbJ",
        "Skew",
        "pmax",
        "HLips",
        "HLi",
        "az",
        "phib",
        "phit",
        "db",
        "dt",
        "lwr",
        "ff",
    ]
    l = [d[p][step] if p in d else None for p in params]
    return l


def plant_parameter(
    surface,
    height,
    pseudostem_height,
    nb_phy,
    nb_young_phy,
    # leaf area dist function
    lad_skew,
    lad_max,
    # internode length distribution
    internode_dist,  # raison geometrique
    sheath_dist,
    # Phyllotaxy for young and mature phytomer
    phyllotactic_angle,
    phyllotactic_young,
    phyllotactic_mature,
    # insertion angle
    basal_insertion,
    deviation_insertion,
    diam_base,
    diam_top,
    leaf_width_ratio,
    leaf_database,
):
    """ """
    nb_phy = int(nb_phy)
    nb_young_phy = int(nb_young_phy)

    relative_phytomer_num = numpy.linspace(1.0 / nb_phy, 1, nb_phy)

    # compute the leaf surface
    k = lad_skew
    nm = lad_max
    r = relative_phytomer_num

    relative_surface = numpy.exp(-k / nm * (2 * (r - nm) ** 2 + (r - nm) ** 3))
    leaf_area = relative_surface / relative_surface.sum() * surface

    # compute normalised surface from database
    db = leaf_database
    rank_max = max(int(x) for x in db)

    def choose_rank(rank):
        if str(rank) not in db:
            rank = rank + 1 if str(rank + 1) in db else rank - 1
        if str(rank) in db:
            return rank

        while str(rank) not in db and rank < rank_max:
            rank += 1
        if str(rank) not in db:
            rank = rank_max
        return rank

    # lindex = ranks = [choose_rank(round(x) * rank_max) for x in relative_phytomer_num]
    norm_surface = numpy.array([fit2(*db[str(rank)][0])[1] for rank in ranks])

    lengths = numpy.sqrt(leaf_area / norm_surface / leaf_width_ratio)
    widths = lengths * leaf_width_ratio

    # internode length
    q = internode_dist
    offset = pseudostem_height
    n = nb_phy - nb_young_phy
    if q == 1:
        u0 = (height - offset) / n
    else:
        u0 = (height - offset) * (1 - q) / (1 - q ** (n + 1))

    internode_length = [0] * nb_young_phy + [u0 * q**i for i in range(n)]

    # sheath length
    q = sheath_dist
    ny = nb_young_phy
    if q == 1:
        u0 = pseudostem_height / ny
    else:
        u0 = pseudostem_height * (1 - q) / (1 - q ** (ny + 1))

    sheath_length = numpy.concatenate(
        (
            numpy.array([u0 * q**i for i in range(ny)]).cumsum(),
            [pseudostem_height] * (nb_phy - ny),
        )
    )

    # internode and sheath diameters
    diameters = [diam_base] * nb_young_phy + numpy.linspace(
        diam_base, diam_top, n
    ).tolist()

    # setup of dictionaries to build  setAdel's dataframes
    NaN = float("nan")

    # axe table
    phyl = 110.0
    startdate = 0
    tip_lig_delay = 1.6 * phyl
    lig_senescence_delay = 10000  # 3.6*phyl for a simulation with realistic senescence
    senescence_disparition_delay = 200

    axeTags = "plant,axe,nf,end,disp,dimIndex,phenIndex,earIndex,emf1,ligf1,senf1,dispf1".split(
        ","
    )
    axeVals = [
        1,
        0,
        nb_phy,
        NaN,
        NaN,
        1,
        1,
        1,
        startdate,
        startdate + tip_lig_delay,
        startdate + tip_lig_delay + lig_senescence_delay,
        startdate + tip_lig_delay + lig_senescence_delay + senescence_disparition_delay,
    ]
    axeT = dict(list(zip(axeTags, axeVals)))

    # phen table
    tip = numpy.array([-phyl, phyl * nb_phy])
    phenTags = ["index", "nrel", "tip", "col", "ssi", "disp"]
    phenVals = [
        [1] * 2,
        [0, 1],
        tip,
        tip + tip_lig_delay,
        tip + tip_lig_delay + lig_senescence_delay,
        tip + tip_lig_delay + lig_senescence_delay + senescence_disparition_delay,
    ]
    phenT = dict(list(zip(phenTags, phenVals)))

    # dim Table
    dimTags = [
        "index",
        "nrel",
        "Ll",
        "Lw",
        "Gl",
        "Gd",
        "El",
        "Ed",
        "incB",
        "dincB",
        "pAngle",
        "dpAngle",
    ]  # last two columns to be tested for existence in setAdel for freeing of geoLeaf function inputs
    dimVals = [
        [1] * nb_phy,
        r,
        lengths,
        widths,
        sheath_length,
        diameters,
        internode_length,
        diameters,
        [basal_insertion] * nb_phy,
        [deviation_insertion] * nb_phy,
        [phyllotactic_angle] * nb_phy,
        [phyllotactic_young] * nb_young_phy
        + [phyllotactic_mature] * (nb_phy - nb_young_phy),
    ]

    dimT = dict(list(zip(dimTags, dimVals)))

    devT = dict(dimT=dimT, phenT=phenT, axeT=axeT, earT=None, ssisenT=None)

    # conversion en liste R de dataframe pour setAdel
    d = dict((k, dataframe(v)) for k, v in devT.items())
    RdevT = robjects.r["list"](**d)
    # RdevT = devT
    return devT, RdevT


# def setAdel(*args, **kwds):
#   return AdelR.setAdel(*args, **kwds),


def setCanopy(canT, *args, **kwds):
    if isinstance(canT, dict):
        canT = dataframe(canT)
    return AdelR.setCanopy(canT, *args, **kwds)
