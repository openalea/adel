def calculateTCViaPovray(
    plantgl_scene,
    output_directory="",
    thermal_time=0,
    camera_distance=200.0,
    fov=50.0,
    width=4288,
    height=2848,
    domain=((0, 0), (1, 1)),
    azimuth=0.0,
    zenith=0.0,
    camera_type="perspective",
    soil=False,
    command="povray",
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
    # genFilepath
    import os

    pov_file = os.path.join(
        output_directory, "%s%s%s" % ("scene_", thermal_time, ".pov")
    )
    # povray
    from openalea.adel.povray import povray

    povray_image, stand_box_image = povray.povray(
        plantgl_scene,
        pov_file,
        camera_distance,
        fov,
        width,
        height,
        domain,
        azimuth,
        zenith,
        camera_type,
        soil,
        command,
    )
    # count_pixels
    count_pixels_result = os.path.join(
        output_directory, "%s%s" % ("_outputpovray_", ".csv")
    )
    from openalea.adel.povray import post_processing

    count_pixels_result = post_processing.count_pixels(
        povray_image, stand_box_image, colors, count_pixels_result
    )
    return (count_pixels_result,)
