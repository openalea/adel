# -*- python -*-
#
#       adel.povray
#
#       Copyright 2006-2012 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#                       Bruno Andrieu <bruno.andrieu@grignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################

import pandas
import numpy as np
from PIL import Image
from pathlib import Path as path

WHITE = (255, 255, 255)


def count_pixels(
    scene_path="",
    scene_box_path="",
    rgb_colors=[],
    result_path="",
    color_labels=[],
    normalize=False,
    absolute_tolerance=1,
    normalized_scene_path=None,
):
    """Counts the pixels of each color, after normalizing the colors or not.
    :Parameters:
        - `scene_path` : the path of the bmp file which represents the scene
        image. MANDATORY.
        - `scene_box_path` : the path of the bmp file which represents the area
        where the pixels have to be counted. MANDATORY.
        - `rgb_colors` : the classes of colors to consider when counting the pixels.
        Each class is described by 3 values: R, G and B. MANDATORY.
        - `result_path` : the path of csv file where the result is added. MANDATORY.
        - `color_labels` : the labels associated to rgb_colors. If empty
        (default), auto-generated labels are: ['class 1', 'class 2', ...]. The
        labels must be given in the same order as rgb_colors, and contains either
        the same number of elements or no element.
        - `normalize` : if True, normalize to 255 the scene and the colors
        before counting the pixels. The colors of the normalized scene which don't
        correspond to the normalized colors given in input are then set to
        (255, 255, 255). Default is False.
        - `absolute_tolerance` : the absolute tolerance used when counting the
        number of pixels of each color. Default is 1.
        - `normalized_scene_path` : the path of the bmp file which represents
        the image of the normalized scene. If empty (default), do not save the
        image of the normalized scene.
    :Types:
        - `scene_path` : str
        - `scene_box_path` : str
        - `rgb_colors` : list of list of int
        - `result_path` : str
        - `color_labels` : list of str
        - `normalize` : bool
        - `absolute_tolerance` : int
        - `normalized_scene_path` : str

    :return: the path of the csv file where the result is added.
        The first column contains the base name of the bmp file which
        represents the scene image (i.e. scene_path). Other
        columns are the list of colors to consider. The last column contains
        the number of pixels of the scene box to consider.
    :rtype: str
    """
    assert (
        scene_path != ""
        and scene_box_path != ""
        and len(rgb_colors) != 0
        and result_path != ""
    )
    if len(color_labels) != 0:
        assert len(color_labels) == len(rgb_colors)

    rgb_colors = [tuple(color) for color in rgb_colors]
    rgb_colors = list(set(rgb_colors))
    color_labels = list(set(color_labels))

    def normalize_to_255(rgb):
        return rgb / rgb.max() * 255.0

    if normalize:
        colors_float = tuple(np.array(rgb_colors).astype(float))
        normalized_colors = [
            tuple(normalize_to_255(color).astype(int)) for color in colors_float
        ]
        colors_mapping = dict(list(zip(normalized_colors, rgb_colors)))
        assert len(normalized_colors) == len(set(normalized_colors))

    scene_path = path(scene_path)
    scene_box_path = path(scene_box_path)
    result_path = path(result_path)

    scene_image = Image.open(scene_path)
    scene_box_image = Image.open(scene_box_path)

    scene_array = scene_image.load()
    scene_box_array = scene_box_image.load()
    width, height = scene_image.size

    normalized_scene_image = Image.new(scene_image.mode, scene_image.size)
    normalized_scene_image_array = normalized_scene_image.load()

    pixels_init_number = np.zeros_like(list(range(len(rgb_colors))))
    pixels_numbers = dict(list(zip(rgb_colors, pixels_init_number)))

    scene_box_pixels_number = 0

    for x in range(width):
        for y in range(height):
            if scene_box_array[x, y] == WHITE:
                scene_box_pixels_number += 1
                scene_color = scene_array[x, y]
                if scene_color == (0, 0, 0):
                    continue
                if normalize:
                    scene_color_float = np.array(scene_color).astype(float)
                    normalized_scene_color = normalize_to_255(scene_color_float)
                    normalized_scene_color = normalized_scene_color.astype(int)
                    normalized_scene_color = tuple(normalized_scene_color)
                    normalized_scene_image_array[x, y] = normalized_scene_color
                    found = False
                    for normalized_color in normalized_colors:
                        if np.allclose(
                            np.array(normalized_color),
                            np.array(normalized_scene_color),
                            0,
                            absolute_tolerance,
                        ):
                            pixels_numbers[colors_mapping[normalized_color]] += 1
                            found = True
                            break
                    if not found:
                        normalized_scene_image_array[x, y] = (255, 255, 255)
                else:
                    if scene_color in rgb_colors:
                        pixels_numbers[scene_color] += 1

    if normalized_scene_path is not None:
        normalized_scene_image.save(path(normalized_scene_path), "BMP")

    if result_path.exists():
        result_df = pandas.read_csv(result_path)
    else:
        columns = ["Filename"]
        if len(color_labels) == 0:
            color_labels = ["class %d" % i for i in range(len(rgb_colors))]
        for rgb, label in zip(rgb_colors, color_labels):
            columns.append("%s - %s" % (rgb, label))
        columns.append("scene_box")
        result_df = pandas.DataFrame(columns=columns)

    data = [
        [scene_path.name]
        + [pixels_numbers[color] for color in rgb_colors]
        + [scene_box_pixels_number]
    ]
    new_index = [result_df.index.size]
    new_df = pandas.DataFrame(data, index=new_index, columns=result_df.columns)

    result_df = pandas.concat([result_df, new_df])
    result_df.to_csv(result_path, na_rep="NA", index=False)

    return str(result_path)
