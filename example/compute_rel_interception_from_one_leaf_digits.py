"""Tutorial: Reconstruct a canopy from digitized data and compute the star
per metamer"""

import pandas
import numpy
from os.path import join as pj
import time
import imp
import os

tps1 = time.clock()

FINAL_NUMBER_OF_PLANTS = 1
INTER_ROW = 0.175
NUMBER_OF_DIRECTIONS = 1  # can be either 1, 16 or 46
NUMBER_OF_LONGITUDINAL_DIVISIONS = 1
DATA_DIR = "data/one_flat_leaf"
DATA_IN_DIR = pj(DATA_DIR, "in")
DATA_OUT_DIR = pj(DATA_DIR, "out")

mean_star_per_phytomer_filepath = pj(DATA_OUT_DIR, "mean_star_per_phytomer.csv")

(
    dimT_filepath,
    earT_filepath,
    lamina2D_filepath,
    laminaCurv_filepath,
    params_filepath,
) = list(
    map(
        pj,
        [DATA_IN_DIR] * 5,
        ["dimT.csv", "earT.csv", "lamina2D.RData", "laminaCurv.RData", "params.py"],
    )
)

params = imp.load_source("params", params_filepath)

(dimT, earT) = list(map(pandas.read_csv, (dimT_filepath, earT_filepath)))

earT_grouped = earT.groupby(["Code", "id_plt"])

for Linc in (0, 0.5, 1):
    angle_of_incidence = Linc * 90
    print("angle_of_incidence", angle_of_incidence)

    metamers_parameters = []

    id_plt_global = 1

    for (code, id_plt), dimT_group in dimT.groupby(["Code", "id_plt"]):
        earT_group = earT_grouped.get_group((code, id_plt))
        earT_series = earT_group.ix[earT_group.first_valid_index()]

        L_ear = earT_group["L_ear"][earT_group.first_valid_index()]
        L_peduncle = earT_group["L_peduncle"][earT_group.first_valid_index()]
        W_peduncle = earT_group["W_peduncle"][earT_group.first_valid_index()]
        A_ear = earT_group["A_ear"][earT_group.first_valid_index()]
        L_spike = earT_group["L_spike"][earT_group.first_valid_index()]

        index_phytomer_ntop = [1]
        metamers_parameters_group = pandas.DataFrame(
            columns=dimT.columns, index=list(range(len(index_phytomer_ntop)))
        )
        metamers_parameters_group["id_plt"] = id_plt_global
        metamers_parameters_group["index_phytomer_ntop"] = index_phytomer_ntop
        metamers_parameters_group["L_blade"] = dimT_group["L_blade"].tolist()
        metamers_parameters_group["W_blade"] = dimT_group["W_blade"].tolist()
        metamers_parameters_group["L_sheath"] = dimT_group["L_sheath"].tolist()
        ear_widths = [W_peduncle, A_ear / L_ear, A_ear / L_ear]
        metamers_parameters_group["W_sheath"] = dimT_group["W_sheath"].tolist()
        metamers_parameters_group["L_internode"] = dimT_group["L_internode"].tolist()
        metamers_parameters_group["W_internode"] = dimT_group["W_internode"].tolist()

        # compute insertion height of leaves
        metamers_parameters_group["h_insertion"] = (
            metamers_parameters_group["L_internode"].cumsum()
            + metamers_parameters_group["L_sheath"]
        )
        # keep only upper elements: metamers (4, 3, 2, 1) and ears (0, -1, -2)
        metamers_parameters_group = metamers_parameters_group[
            metamers_parameters_group.index_phytomer_ntop <= 4
        ]
        # create pseudo-stem elements (spacer between leaves)
        metamers_parameters_group["pstem"] = numpy.diff(
            [0] + metamers_parameters_group["h_insertion"].tolist()
        )

        # keep only the useful columns
        metamers_parameters_group = metamers_parameters_group[
            ["id_plt", "index_phytomer_ntop", "L_blade", "W_blade", "pstem", "W_sheath"]
        ]
        # rename the columns using the current convention of adel
        metamers_parameters_group.rename(
            columns={
                "id_plt": "plant",
                "index_phytomer_ntop": "ntop",
                "L_blade": "Ll",
                "W_blade": "Lw_shape",
                "pstem": "El",
                "W_sheath": "Ed",
            },
            inplace=True,
        )

        # add mandatory topological info
        metamers_parameters_group["axe_id"] = (
            "MS"  # Note: to adapt if tillers are provided
        )
        metamers_parameters_group["ms_insertion"] = (
            0  # Note: to adapt if tillers are provided
        )
        # order numphy from base to top
        metamers_parameters_group["numphy"] = list(
            range(1, len(metamers_parameters_group) + 1)
        )
        # add missing mandatory data  (does like adel)
        metamers_parameters_group["Laz"] = [
            180 for i in range(len(metamers_parameters_group))
        ]  # leaf azimuth
        metamers_parameters_group["LcType"] = metamers_parameters_group[
            "ntop"
        ]  # selector for first level in leaf db
        metamers_parameters_group["Lc_index"] = (
            1  # selector for second level in leaf_db (ranging 1:max_nb_leaf_per_level)
        )
        # fill other columns
        metamers_parameters_group["Lv"] = metamers_parameters_group["Ll"]
        metamers_parameters_group["Lr"] = 0
        metamers_parameters_group["Lsen"] = 0
        metamers_parameters_group["L_shape"] = metamers_parameters_group["Ll"]
        metamers_parameters_group["Linc"] = Linc
        metamers_parameters_group["Gv"] = 0
        metamers_parameters_group["Gd"] = 0
        metamers_parameters_group["Gsen"] = 0
        metamers_parameters_group["Ginc"] = 0
        metamers_parameters_group["Ev"] = metamers_parameters_group["El"]
        metamers_parameters_group["Esen"] = 0
        metamers_parameters_group["Einc"] = 0  # Note: to adapt if azimuth is provided

        metamers_parameters.append(metamers_parameters_group)

        id_plt_global += 1

    number_of_available_plants = len(metamers_parameters)

    if number_of_available_plants < FINAL_NUMBER_OF_PLANTS:
        # sample the available plants to have FINAL_NUMBER_OF_PLANTS in the constructed canopy
        number_of_plants_to_sample = FINAL_NUMBER_OF_PLANTS - number_of_available_plants
        metamers_parameters_sample = numpy.random.choice(
            metamers_parameters, number_of_plants_to_sample
        )
        metamers_parameters.extend(metamers_parameters_sample)

    metamers_parameters = pandas.concat(metamers_parameters, ignore_index=True)

    metamers_parameters_filepath = pj(
        DATA_OUT_DIR, "metamers_parameters" + str(angle_of_incidence) + ".csv"
    )
    print("Save Adel parameters to ", metamers_parameters_filepath)
    metamers_parameters.to_csv(metamers_parameters_filepath, na_rep="NA", index=False)

    def agronomicplot(
        nplants=None,
        length=None,
        width=None,
        sowing_density=None,
        plant_density=None,
        inter_row=None,
        noise=0,
        convunit=100,
        center_scene=True,
    ):
        """Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
        nplants is the number of plants of the micro-plot (mandatory if length or width is None)
        length (m) is plot dimension along row direction (mandatory if nplants is None)
        width (m) is plot dimension perpendicular to row direction (mandatory if nplants is None)
        sowing density is the density of seeds sawn (mandatory)
        plant_density is the density of plants that are present (after loss due to bad emergence, early death...) (mandatory)
        inter_row (m) is for the  distance between rows (mandatory)
        noise (%), indicates the precision of the sowing for the inter plant spacing
        convunit is the conversion factor from meter to scene unit
        center_scene allows to center the position around origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain

        Rows are parallel to x-axis
        Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
        """
        from operator import itemgetter
        import random
        from openalea.adel.stand.stand import regular

        inter_plant = 1.0 / inter_row / sowing_density

        if nplants is None:  # length and width are NOT None
            nrow = max(1, int(float(width) / inter_row))
            plant_per_row = max(1, int(float(length) / inter_plant))
            nplants = nrow * plant_per_row

        else:  # nplants is NOT None ; length and width are None
            import math

            nrow = max(1, int(math.sqrt(nplants)))

        dx = inter_plant * convunit
        dy = inter_row * convunit
        positions, domain = regular(nplants, nrow, dx, dy)
        n_emerged = int(nplants * plant_density / sowing_density)
        positions = random.sample(positions, n_emerged)
        # sorting by ranks
        positions = sorted(positions, key=itemgetter(1, 0))
        domain_area = (
            abs(domain[1][0] - domain[0][0])
            / convunit
            * abs(domain[1][1] - domain[0][1])
            / convunit
        )
        if center_scene:
            xc = float(domain[1][0] + domain[0][0]) / 2
            yc = float(domain[1][1] + domain[0][1]) / 2
            positions = [(x - xc, y - yc, z) for x, y, z in positions]
            domain = (
                (domain[0][0] - xc, domain[0][1] - yc),
                (domain[1][0] - xc, domain[1][1] - yc),
            )

        return n_emerged, positions, domain, domain_area, 1.0 / convunit

    # create a canopy with the same number of plants as in the canopy table (mandatory)
    nplants, positions, domain, domain_area, convUnit = agronomicplot(
        nplants=FINAL_NUMBER_OF_PLANTS,
        sowing_density=params.density,
        plant_density=params.density,
        inter_row=INTER_ROW,
    )

    # create mtg

    from alinea.adel.newmtg import mtg_factory, adel_metamer, adel_label
    from alinea.adel.mtg_interpreter import mtg_interpreter, plot3d
    from openalea.plantgl.all import Viewer

    ## geometry
    # extract_leaves
    from alinea.adel.wheat import extract_wheat

    db = extract_wheat.extract_leaf_info(
        pj(DATA_IN_DIR, "laminaCurv.RData"), pj(DATA_IN_DIR, "lamina2D.RData")
    )
    # fit leaves
    from alinea.adel.fit import fit

    db, discard = fit.fit_leaves(db, NUMBER_OF_LONGITUDINAL_DIVISIONS)

    # set the rotation angle of the plant on itself. Permits to randomize the orientation of each leaf.
    stand = list(zip(positions, [0 for i in range(len(positions))]))
    g = mtg_factory(
        metamers_parameters.to_dict("list"), adel_metamer, leaf_db=db, stand=stand
    )
    # add geometry
    g = mtg_interpreter(g)

    ## plot (uncomment the next 2 lines to plot the 3D model)
    scene = plot3d(g)
    Viewer.display(scene)

    # call caribu
    from alinea.caribu.caribu_star import caribu_star

    geom = g.property("geometry")
    star, exposed_area = caribu_star(
        geom, directions=NUMBER_OF_DIRECTIONS, domain=domain, convUnit=convUnit
    )  # cf caribu_star doc for output interpretation
    caribu_out = pandas.DataFrame(
        [(adel_label(g, vid), star[vid], exposed_area[vid]) for vid in star],
        columns=["adel_label", "star", "exposed_area"],
    )

    caribu_out_filepath = pj(
        DATA_OUT_DIR, "caribu_out" + str(angle_of_incidence) + ".csv"
    )
    print("Save caribu output to", caribu_out_filepath)
    caribu_out.to_csv(caribu_out_filepath, na_rep="NA", index=False)

    # Calculate the star for each class of metamers (metamer1, metamer2, metamer3, etc).
    # Star is "ratio viewed area / area". It takes into consideration the angles of
    # incidence ('directions') and the hidden fractions of the leaves
    # Extract the phytomer number for each element in res['adel_label']:
    last_part_of_adel_label = numpy.char.partition(
        caribu_out["adel_label"].values.astype(str), "metamer"
    )[:, 2]
    metamers = numpy.char.partition(last_part_of_adel_label, "_")[:, 0]
    # create a column 'numphy' in the dataframe caribu_out
    caribu_out["numphy"] = metamers
    caribu_out["angle_of_incidence"] = angle_of_incidence

    mean_star_per_phytomer = caribu_out[["angle_of_incidence", "star"]]

    if os.path.exists(mean_star_per_phytomer_filepath):
        mean_star_per_phytomer_pre = pandas.read_csv(mean_star_per_phytomer_filepath)
        mean_star_per_phytomer = pandas.concat(
            [mean_star_per_phytomer_pre, mean_star_per_phytomer], ignore_index=True
        )

    print("Save mean star per phytomer to", mean_star_per_phytomer_filepath)
    mean_star_per_phytomer.to_csv(
        mean_star_per_phytomer_filepath, na_rep="NA", index=False
    )

tps2 = time.clock()
print("Execution time: ", tps2 - tps1)
