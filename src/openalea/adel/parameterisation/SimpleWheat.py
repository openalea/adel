import csv as csv
from math import ceil, floor
from numpy import array, arange
from copy import deepcopy


class dim_pattern:
    """
    Class which implements patterns for predicting dimension, Aera_contrib and so on from data.
    """

    def __init__(self, dimT, index=111):
        """
        Initialisation of dimension pattern, just need the equivalent of a dimT
        dimT:
        index: column which identify axis in dimT
        """
        # Opening dimT file
        csvfile = open(dimT, "rb")
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        # Initialisation of the usefull information
        self.decphy = [0, 3, 4, 5, 6, 7]  # decalage des phytomeres pour les talles
        ndict = {
            "relative_phytomer": [],
            "surface": [],
            "bladelength": [],
            "sheathlength": [],
            "internodelength": [],
            "stemdiameter": [],
        }
        self.emptydict = ndict
        self.dictMS = deepcopy(ndict)
        # dimT construction CAREFULL NO HEADER
        next(reader)  # Getting the header out
        row = next(reader)  # initialisation
        while int(row[0]) is not index:
            row = next(reader)

        while int(row[0]) is index:
            self.dictMS["relative_phytomer"].append(float(row[1]))
            self.dictMS["surface"].append(float(row[2]) * float(row[3]) * float(row[9]))
            self.dictMS["bladelength"].append(float(row[2]))
            self.dictMS["sheathlength"].append(float(row[4]))
            self.dictMS["internodelength"].append(float(row[6]))
            self.dictMS["stemdiameter"].append(float(row[5]))
            row = next(reader)
        self.numphy = len(
            self.dictMS["stemdiameter"]
        )  # number of phytomer for mainstem

    def predict_tiller(self, tiller_position):
        """
        Returns datadict for tiller_position: 0 is main stem, 5 is the absolute last tiller
        """
        tiller_position = int(tiller_position)
        if tiller_position == 0:
            return self.dictMS
        else:
            dict = deepcopy(
                self.emptydict
            )  # Take the empty dictionnary having the same ontology then dictMS
            for i in range(self.decphy[tiller_position], self.numphy):
                dict["relative_phytomer"].append(self.dictMS["relative_phytomer"][i])
                dict["surface"].append(self.dictMS["surface"][i])
                dict["bladelength"].append(self.dictMS["bladelength"][i])
                dict["sheathlength"].append(self.dictMS["sheathlength"][i])
                dict["internodelength"].append(self.dictMS["internodelength"][i])
                dict["stemdiameter"].append(self.dictMS["stemdiameter"][i])
            return dict

    def total_area(self, dict):
        """
        Returns total aera of a tiller which is represented by a dict
        """
        return sum(dict["surface"])

    def PlantArea_contrib(self, nb_tiller):
        """
        Contribution of each tiller for the total area
        """
        surface = []
        for i in arange(nb_tiller + 1):  # We include the main stem
            surface.append(self.total_area(self.predict_tiller(i)))
        Stot = sum(array(surface))
        return array(surface) / Stot
        # out = list_poid_0 a nb tiller


soisson = dim_pattern("../data/dimTSoissons.csv")


def Wheatpop(
    LAI=3.3,
    axe_density=1001,
    plant_density=300.0,
    nout=10,
    nb_phy=11,
    dim_pattern=soisson,  # Qui contient PlantArea_contrib (% de plant area porte par les axes)
):
    """
    Generate a population of nout different axe with appropriate tiller number and surface per plant
    To be plug to "agronomic plot" and "SimpleWheat"

    TO DO :
    - > donner un stade ET un nombre de feuille verte contribuant au LAI sur le bm => on sait ou se poisitionner dans le pattern
    - > donner un LAI qui represente la surface des feuilles presente (vert + sene) et un % senescent global sur ce LAI
    - > + et la hauteur

    LAI
    axe_density
    plant_density
    nout is number of plant
    nb_phy: nb_phytomer
    """
    # How many tiller for each main stem
    target_surface = LAI / plant_density  # target_surface is in m2
    mean_tiller_per_plant = axe_density / plant_density - 1  # Main stem aren't tillers
    tiller_min = floor(mean_tiller_per_plant)

    nb_plant_max = floor(nout * (mean_tiller_per_plant - tiller_min))
    nb_plant_min = nout - nb_plant_max
    # More tillers, more surface => correction
    AreaContrib_tiller_max = dim_pattern.PlantArea_contrib(
        tiller_min + 1
    )  # tiller_max = tiller_min +1
    weight_plant_tiller_min = 1  # reference for min nb of tiller
    weight_plant_tiller_max = 1 + AreaContrib_tiller_max[tiller_min + 1] / (
        1 - AreaContrib_tiller_max[tiller_min + 1]
    )  # Area of last tiller
    total_weight = (
        weight_plant_tiller_min * nb_plant_min + weight_plant_tiller_max * nb_plant_max
    )

    Plant_min = array(
        [
            nb_plant_min,
            tiller_min + 1,
            nb_phy,
            target_surface * weight_plant_tiller_min / total_weight,
        ]
    )
    Plant_max = array(
        [
            nb_plant_max,
            tiller_min + 2,
            nb_phy,
            target_surface * weight_plant_tiller_max / total_weight,
        ]
    )
    # nb de plante, nb_axe, nb_phyt, surface_plante
    return array([Plant_min, Plant_max])


def simpleWheat_param(
    total_area=10000,  # aussi dans wheatpop
    total_height=200,  # Nouveau!!!
    nb_phy=11,  # entree de Wheatpop
    # leaf area dist function
    dim_pattern="soisson",  # function des distributions de longueurs d'organes et du passage brin-maitre tiller(gaine entre-noeud et limbe) (3functions)
    geom_pattern="soisson",  # BD des XYSR de nervures mesurees
    # dred: maximal distance between maintsteem and tiller at flowering
    dred=7,  # Nouveau
    diam_base=3,  # Nouveau
    diam_top=1,  # Nouveau
    nbtalles=3,  # topo: info dans le mtg (mtg.verticies_at_scale(1))
    # No azimuth param
):
    """
    generate simpleMais dictionary from parameter for wheat-like plants

    """

    # creation brin maitre

    nb_phy = int(nb_phy)
    nb_young_phy = int(nb_young_phy)

    # compute the leaf surface
    leaf_area = np.array(bell_shaped_dist(total_area, nb_phy, lad_rmax, lad_skew))
    # derive corresponding length and widths
    lengths = np.sqrt(leaf_area / np.array(shape_factors) / leaf_width_ratio)
    widths = lengths * leaf_width_ratio

    # distances between leaves
    pseudostem = geometric_dist(pseudostem_height, nb_young_phy, pseudostem_dist)
    stem = geometric_dist(
        total_height - pseudostem_height, nb_phy - nb_young_phy, stem_dist
    )
    distances = pseudostem + stem
    # sheath and internode true length
    pseudostem_sheath = (np.array(pseudostem)).cumsum()
    pseudostem_internode = [0] * nb_young_phy
    stem_sheath = np.array(stem)
    stem_internode = [0] + (stem_sheath[1:]).tolist()
    sheaths = pseudostem_sheath.tolist() + stem_sheath.tolist()
    internodes = pseudostem_internode + stem_internode
    # stem diameters
    diameters = [diam_base] * nb_young_phy + np.linspace(
        diam_base, diam_top, nb_phy - nb_young_phy
    ).tolist()
    # inclinations
    inclinations = np.linspace(basal_insertion, top_insertion, nb_phy)
    # azimuths Attention R
    # Azimuth = (n-4 n) 180 + 60 * randn(1), (1 n-4)   180 + 20 * randn(1))

    axe = [0] * nbphy
    Einc = axe
    # creation des axes

    # Alexis: initialisations des talles
    t_lengths = talle_expend_length(lengths, nbtalles, nb_phy)
    t_widths = talle_expend_length(widths, nbtalles, nb_phy)
    t_distances = talle_expend_length(distances, nbtalles, nb_phy)
    t_sheaths = talle_expend_length(sheaths, nbtalles, nb_phy)
    t_internodes = talle_expend_length(internodes, nbtalles, nb_phy)

    for i in range(1, nbtalles + 1):
        t_nb_phy = round(nb_phy - (2.5 + 0.5 * i))
        t_distances = distances[(nb_phy - t_nb_phy - 1) :]
        distances = np.concatenate((distances, t_distances))


        inclinations = np.linspace(basal_insertion, top_insertion, t_nb_phy)
        # azimuth : aller chercher regle adel ?
        axe = axe + [i] * t_nb_phy
        #
        t_Einc = [60] + [0] * (t_nb_phy - 1)
    # output dict
    dTags = [
        "plante",
        "axe",
        "phytomere",
        "longFeu",
        "largFeu",
        "inclinaisonFeu",
        "azimuthFeu",
        "longTige",
        "diamTige",
        "longGa",
        "longEn",
        "incEn",
    ]
    dVals = [
        [1] * nb_phy * (nbtalles + 1),
        axe,
        list(range(1, nb_phy + 1)),
        lengths,
        widths,
        inclinations,
        azimuths,
        distances,
        diameters,
        sheaths,
        internodes,
        Einc,
    ]
    dVals = list(map(np.array, dVals))
    dout = dict(list(zip(dTags, dVals)))

    return dout


# Alexis: propage l'info d'une talle a la suivante
def talle_expend_length(param, nbtalles, nb_phy):
    exp = param  # Initialisation, param appartient uniquement au MB
    for i in range(1, nbtalles + 1):
        t_nb_phy = round(nb_phy - (2.5 + 0.5 * i))
        t_param = param[(nb_phy - t_nb_phy - 1) :]
        exp = np.concatenate((exp, t_param))
    return exp
