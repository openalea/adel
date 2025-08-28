def energy_per_organ(
    nplants,
    positions,
    adel_output,
    leaves_db,
    domain,
    classic=True,
    number_of_directions=46,
    convUnit=0.01,
    caribu_output="Ei",
):
    import numpy as np
    import pandas as pd
    from openalea.adel.newmtg import mtg_factory, adel_metamer, adel_label
    from openalea.adel.mtg_interpreter import mtg_interpreter
    from openalea.caribu.caribu_star import run_caribu, diffuse_source

    # rotation angle of the plant on itself.
    azimuths = np.random.random(nplants) * 360  # TODO: use input
    stand = list(zip(positions, azimuths))
    g = mtg_factory(adel_output, adel_metamer, leaf_db=leaves_db, stand=stand)
    # add geometry
    g = mtg_interpreter(g, classic=classic)
    # call caribu
    geom = g.property("geometry")

    sources = diffuse_source(number_of_directions)
    out = run_caribu(sources, geom, domain=domain)
    selected_out = out[caribu_output]
    areas = out["Area"]
    exposed_area = {vid: selected_out[vid] * areas[vid] * convUnit**2 for vid in areas}

    energy = np.array(
        [
            (adel_label(g, vid), selected_out[vid], exposed_area[vid])
            for vid in selected_out
        ]
    )
    adel_labels_array = np.array(energy[:, 0])
    adel_labels_split = np.array(np.char.split(adel_labels_array, "_").tolist())
    plant_labels = np.char.strip(adel_labels_split[:, 0], "plant")
    axes_labels = adel_labels_split[:, 1]
    metamer_labels = np.char.strip(adel_labels_split[:, 2], "metamer")
    organ_labels = adel_labels_split[:, 3]
    elements_labels = adel_labels_split[:, 4]
    star_array = np.array(energy[:, 1])
    exposed_area_array = np.array(energy[:, 2])

    energy_enlarged = np.array(
        [
            plant_labels,
            axes_labels,
            metamer_labels,
            organ_labels,
            elements_labels,
            star_array,
            exposed_area_array,
        ]
    )

    energy_df = pd.DataFrame(
        energy_enlarged.transpose(),
        columns=[
            "plant",
            "axis",
            "metamer",
            "organ",
            "element",
            caribu_output,
            "exposed_area",
        ],
    )
    energy_df[["plant", "metamer"]] = energy_df[["plant", "metamer"]].astype(int)
    energy_df.sort(["plant", "axis", "metamer", "organ", "element"], inplace=True)

    return (energy_df,)


def energy_per_organ_macro(
    nplants,
    positions,
    adel_output,
    leaves_db,
    domain,
    output_directory,
    thermal_time,
    classic=True,
    number_of_directions=46,
    convUnit=0.01,
    caribu_output="Ei",
):
    energy_df = energy_per_organ(
        nplants,
        positions,
        adel_output,
        leaves_db,
        domain,
        classic,
        number_of_directions,
        convUnit,
        caribu_output,
    )[0]
    import os

    energy_path = os.path.join(
        output_directory, "%s%s%s" % ("_energy_", thermal_time, ".csv")
    )
    energy_df.to_csv(energy_path, na_rep="NA", index=False)
    return (energy_path,)
