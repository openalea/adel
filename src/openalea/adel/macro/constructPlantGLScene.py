def constructPlantGLScene(
    stand_MTG,
    property_name="tissue_type",
    colors=[
        [0, 255, 0],
        [255, 0, 0],
        [0, 255, 0],
        [255, 0, 0],
        [0, 255, 0],
        [255, 0, 0],
        [0, 255, 0],
        [255, 0, 0],
        [0, 255, 0],
        [255, 0, 0],
        [0, 255, 0],
        [255, 0, 0],
    ],
):
    # col_item
    from openalea.adel.povray import povray

    lambda_function = povray.col_item(None, colors)
    # apply_property
    from openalea.adel import mtg

    new_values = mtg.apply_property(stand_MTG, property_name, lambda_function)
    # to_plantgl
    plantgl_scene = mtg.to_plantgl(
        stand_MTG, (0, 180, 0), (0, 130, 0), (170, 85, 0), new_values, True
    )[0]
    return (plantgl_scene,)
