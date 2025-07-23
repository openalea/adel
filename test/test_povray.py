import pandas
import numpy as np
from pathlib import Path as path
from openalea.adel.povray import post_processing


def test_count_pixels():
    scene_path =  path(__file__).parent.resolve()  / path("data/test_povray/scene.bmp")
    scene_box_path =  path(__file__).parent.resolve() / path("data/test_povray/scene_box.bmp")

    rgb_colors = [
        [0, 170, 0],
        [170, 170, 0],
        [255, 0, 0],
        [255, 170, 0],
        [255, 201, 5],
        [0, 0, 255],
        [255, 170, 0],
        [0, 0, 255],
        [255, 201, 5],
        [170, 255, 255],
        [170, 255, 255],
    ]

    import tempfile

    povray_results_directory =  path(__file__).parent.resolve() / path(tempfile.mkdtemp(suffix="_povray_results"))
    result_path = povray_results_directory.joinpath("scene_res.csv")
    normalized_scene_path = povray_results_directory.joinpath("normalized_scene.bmp")
    normalize = True
    post_processing.count_pixels(
        scene_path=scene_path,
        scene_box_path=scene_box_path,
        rgb_colors=rgb_colors,
        result_path=result_path,
        normalize=normalize,
        normalized_scene_path=normalized_scene_path,
    )
    expected_result_array = np.array(
        [["scene.bmp", 0, 28526, 375, 2253, 25, 83, 2959, 122786]], dtype=object
    )
    result_array = pandas.read_csv(result_path).values
    np.testing.assert_equal(result_array, expected_result_array)
