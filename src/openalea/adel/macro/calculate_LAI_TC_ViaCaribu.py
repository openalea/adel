def calculate_LAI_TC_ViaCaribu(
    stand_MTG, pattern=None, thermal_time=0, output_directory=""
):
    # light string
    from openalea.caribu import lightString

    directional_vector = lightString.lightString(1.0, 0.0, 46.0)
    # par4.opt
    from openalea.core import pkgmanager

    pm = pkgmanager.PackageManager()
    pm.init(verbose=False)
    from openalea.core.system import systemnodes

    opt_filepath = list(
        systemnodes.get_data("par4.opt", "openalea.caribu.data").values()
    ).pop()
    # read
    from openalea.file import files

    opt_str = files.FileRead()(opt_filepath)
    # CaribuScene
    from openalea.caribu import CaribuScene_nodes

    caribu_scene = CaribuScene_nodes.newCaribuScene(
        stand_MTG, directional_vector, pattern, opt_str
    )[0]
    # addsoil_2
    from openalea.adel.macro import addsoil_2

    caribu_scene = addsoil_2.addsoil_2(caribu_scene, 2.0)
    # Caribu
    caribu_scene, energy = CaribuScene_nodes.runCaribu(
        caribu_scene,
        True,
        {"SphereDiameter": 0.5, "Nz": 5, "Zmax": 2, "keepFF": False},
        True,
    )
    # LIE
    from openalea.caribu import selectOutput

    canestra_output_1 = selectOutput.selectOutput(energy, "Opt")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Einc")[0]
    from openalea.caribu import filterby

    filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x == 0.0
    )[0]
    total_soil = sum(filtered)
    total_incident = CaribuScene_nodes.getIncidentEnergy(caribu_scene)[2]
    efficience = (total_incident - total_soil) / float(total_incident)
    # LAIg
    canestra_output_1 = selectOutput.selectOutput(energy, "Opt")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Area")[0]
    soil_filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x == 0.0
    )[0]
    soil_area = sum(soil_filtered)

    canestra_output_1 = selectOutput.selectOutput(energy, "Opt")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Area")[0]
    canopy_triangles_filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x > 0.0
    )[0]
    canopy_triangles_area = sum(canopy_triangles_filtered)

    canestra_output_1 = selectOutput.selectOutput(energy, "Opt")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Area")[0]
    green_triangles_filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x == 1.0
    )[0]
    green_triangles_area = sum(green_triangles_filtered)

    canestra_output_1 = selectOutput.selectOutput(energy, "Opak")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Area")[0]
    opak_leaves_filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x > 0.0
    )[0]
    opak_leaves_area = sum(opak_leaves_filtered)

    canestra_output_1 = selectOutput.selectOutput(energy, "Opak")[0]
    canestra_output_2 = selectOutput.selectOutput(energy, "Opt")[0]
    opak_opt_filtered = filterby.filterby(
        canestra_output_1, canestra_output_2, lambda x: x > 0.0
    )[0]

    LAI = opak_leaves_area / float(soil_area)
    from openalea.adel.macro import AND

    LAIg = sum(AND.AND(opak_leaves_filtered, opak_opt_filtered)[0]) / float(soil_area)
    CAI = (canopy_triangles_area / float(soil_area) - LAI) / 2 + LAI
    CAIg = (green_triangles_area / float(soil_area) - LAIg) / 2 + LAIg

    # genFilepath
    import os

    LAI_PAI_path = os.path.join(output_directory, "%s%s" % ("_caribuoutput_", ".csv"))

    # Write_LAI_PAI
    from openalea.adel.macro import write_LAI_PAI

    LAI_PAI_path = write_LAI_PAI.Write_LAI_PAI(
        thermal_time, efficience, CAI, CAIg, LAI, LAIg, LAI_PAI_path
    )

    return (LAI_PAI_path,)
