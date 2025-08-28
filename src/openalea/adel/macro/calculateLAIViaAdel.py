def calculateLAIViaAdel(
    canopy_table,
    Lstrings,
    lamina2D_filepath="",
    thermal_time=0,
    output_directory="",
    number_of_plants=0,
    domain_area=0.0,
):
    # canL2canS
    from openalea.adel import AdelR

    canopy_table = AdelR.canL2canS(canopy_table, lamina2D_filepath, 0.5)

    # update
    Lstrings.update(canopy_table)

    # Write_Table
    from openalea.adel.wheat import Write_Table
    import os

    canl2canS_filepath = os.path.join(
        output_directory, "%s%s%s" % ("_canl2canS_", thermal_time, ".csv")
    )
    adel_output_filepath = Write_Table.Write_Table(
        Lstrings, canl2canS_filepath, sep=","
    )[0]

    # post_processing
    global_postprocessing_filepath = os.path.join(
        output_directory, "%s%s" % ("_postprocessing_global", ".csv")
    )
    peraxis_postprocessing_filepath = os.path.join(
        output_directory, "%s%s" % ("_postprocessing_peraxis_", ".csv")
    )
    from openalea.adel.stand import stand

    (
        global_postprocessing_filepath,
        peraxis_postprocessing_filepath,
        intermediate_postprocessing_filepath,
    ) = stand.post_processing(
        adel_output_filepath,
        number_of_plants,
        domain_area,
        global_postprocessing_filepath,
        peraxis_postprocessing_filepath,
    )

    return (
        global_postprocessing_filepath,
        peraxis_postprocessing_filepath,
        intermediate_postprocessing_filepath,
    )
