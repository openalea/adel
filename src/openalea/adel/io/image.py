from openalea.plantgl.all import Viewer
import os.path


def save_image(scene, image_name="%s/img%04d.%s", directory=".", index=0, ext="png"):
    """
    Save an image of a scene in a specific directory

    Parameters
    ----------

        - scene: a PlantGL scene
        - image_name: a string template
            The format of the string is dir/img5.png
        - directory (optional: ".") the directory where the images are written
        - index: the index of the image
        - ext : the image format

    Example
    -------

        - Movie:
            convert *.png movie.mpeg
            convert *.png movie.gif
            mencoder "mf://*.png" -mf type=png:fps=25 -ovc lavc -o output.avi
            mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts vcodec=msmpeg4v2:vhq "mf://*.png" -mf type=png:fps=18 -of avi  -o output.avi

    """

    if not image_name:
        image_name = "{directory}/img{index:0>4d}.{ext}"
    filename = image_name.format(directory=directory, index=index, ext=ext)
    Viewer.frameGL.saveImage(filename)
    return (scene,)


# def movie(directory, input_file, output_file, fps=25, encoder='mencoder'):
