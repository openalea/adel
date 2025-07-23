# -*- python -*-
#
#       adel.povray
#
#       Copyright 2006-2014 INRIA - CIRAD - INRA
#
#       File author(s): Christian Fournier <christian.fournier@supagro.inra.fr>
#                       Christophe Pradal <christophe.pradal@inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
import os
import tempfile
import platform
from pathlib import Path as path
import openalea.plantgl.all as pgl
from openalea.adel.postprocessing import domain3D, stand_box


class PovRayError(Exception):
    pass


class PovRay:
    """A class interface to povray Raytracer"""

    def __init__(
        self,
        camera={
            "type": "perspective",
            "distance": 1.0,
            "fov": 45.0,
            "xc": 0.0,
            "yc": 0.0,
            "azimuth": 0,
            "zenith": 0.0,
        },
        image_width=320,
        image_height=280,
        background=(0, 0, 0),
        light_position=(10, 0, 10),
        light_color=(1, 1, 1),
        working_dir=None,
    ):
        """Setup a Povray instance

        :Parameters:
            - camera: a dict of parameters for positioning the camera
                - distance: distance from the position of the camera to the look_at point
                - fov : angle corresponding to the width of the image
                - xc, yc : coordinates of the center of the scene
                - azimuth: angle in degree around the vertical axis. az=0. is equialent to have the width of the image align to X direction of the scene
                - zenith: angle between the view direction and the vertical
            - image_width: width of the final image in pixel
            - image_height: height of the final image in pixel
            - working_dir: A directory for storing povary files. If None, (default),  a temporary directory will be created and removed onece the instance is deleted

        """

        try:
            if working_dir is not None:
                self.wdir = path(path(working_dir).abspath())
                if not self.wdir.exists():
                    self.wdir.mkdir()
                self.cleanup_wdir = False
            else:
                # build a temporary directory
                self.wdir = path(tempfile.mkdtemp())
                self.cleanup_wdir = True
        except:
            raise PovRayError(
                "PovRay can't create its working directory : check for read/write permission or security level"
            )

        if platform.system() == "Windows":
            self.cmdline = "pvengine +FN +I%s +H%d +W%d -d /exit"
        else:
            self.cmdline = "povray +FN +I%s +H%d +W%d"

        self.image_width = image_width
        self.image_height = image_height
        self.camera = camera
        self.light_position = light_position
        self.light_color = light_color
        self.background = background

        self.soil = False
        self.rendered_image_path = None
        self.tesselator = pgl.Tesselator()

        # to be removed once cv_camera has been properly integrated
        self.user_camera = None

    def __del__(self):
        if self.cleanup_wdir:
            if self.wdir.exists():
                # print 'Remove tempfile %s'%self.wdir
                self.wdir.rmtree()
        else:
            print("Povray.__del__ called, but working directory kept: %s" % self.wdir)

    def add_soil(self, domain):
        self.domain = domain
        self.soil = True

    def camera_string(
        self,
        type="perspective",
        distance=1.0,
        fov=45.0,
        xc=0.0,
        yc=0.0,
        azimuth=0,
        zenith=0.0,
    ):
        """String representation of the camera"""

        pov_camera = """
camera {{
    {camera}
    location <{tx:.2f},{ty:.2f},{tz:.2f}>
    direction <0,0,-1>
    right <{width:d}/{height:d},0,0>
    angle {fov:n}

    rotate <{zenith:.2f},0,0>
    rotate <0,0,{azimuth:.2f}>
}}

    """.format(
            camera=type,
            tz=distance,
            fov=fov,
            azimuth=azimuth,
            zenith=zenith,
            tx=xc,
            ty=yc,
            width=self.image_width,
            height=self.image_height,
        )

        # hack
        if self.user_camera is not None:
            pov_camera = self.user_camera

        return pov_camera

    def set_user_camera(self, pov_camera):
        self.user_camera = pov_camera

    def soil_string(self):
        (x1, y1), (x2, y2) = self.domain
        s = """
#declare T1=
    texture {{
      pigment {{ color red 1.0 green 0.75 blue 0.33 }}
      normal {{bumps 0.4  scale 0.007}} 
}}

#declare T2=
 texture{{
   pigment{{color red 1.0 green 0.75 blue 0.33}}
   normal {{granite 0.6 scale 0.1 }}
   finish {{phong 0.8 phong_size 200}} 
}}

box {{ <{x1}, {y1},  -0.1>,  
    < {x2}, {y2}, 0>   
    texture {{T2}} }}
"""
        return s.format(x1=x1, y1=y1, x2=x2, y2=y2)

    def render(self, scene, name="scene.pov"):
        """Render the scene"""

        old_dir = os.path.abspath(os.getcwd())
        self.rendered_image_path = None

        try:  # catch error to return to original directory
            os.chdir(self.wdir)
            f = self.wdir / name
            namebase = f.namebase

            mesh_fn = self.wdir / (namebase + "_mesh.pov")
            pov = pgl.PovFilePrinter(str(mesh_fn), self.tesselator)
            scene.apply(pov)

            image_name = self.wdir / (namebase + ".png")
            if image_name.exists():
                image_name.remove()

            povf = f.open(mode="w")
            povf.write("/" * 80 + "\n")
            povf.write('#include "colors.inc"' + "\n\n")
            povf.write(self.camera_string(**self.camera) + "\n")
            x, y, z = self.light_position
            r, g, b = self.light_color
            povf.write(
                "light_source {{  <{tx},{ty},{tz}> color rgb <{r},{g},{b}>}}\n\n".format(
                    tx=x, ty=y, tz=z, r=r, g=g, b=b
                )
            )
            r, g, b = self.background
            povf.write(
                "background {{ color rgb < {cr},{cg},{cb}> }}\n\n".format(
                    cr=r, cg=g, cb=b
                )
            )
            if self.soil:
                povf.write(self.soil_string())
            povf.write('#include "{mesh_fn}"\n'.format(mesh_fn=mesh_fn))
            povf.close()

            cmd = self.cmdline % (str(f), self.image_height, self.image_width)
            os.system(cmd)

            if image_name.exists():
                self.rendered_image_path = str(self.wdir / namebase + ".png")
            else:
                raise PovRayError(
                    "Image not created. Check that povray is installed and set in the path (path to the pvengine programm)."
                )
        finally:
            os.chdir(old_dir)

    def get_image(self, reader, **args):
        """return the last rendered image object

        call reader(image_path, **args) to read the image
        """

        im = None
        if self.rendered_image_path is not None:
            im = reader(self.rendered_image_path, **args)
        return im


# deprecated functions (here for backward compatibility)


color_list = [
    (0, 0, 0),
    (255, 0, 0),  # 1 = Green lamina
    (0, 255, 0),  # 2 = Senescent Lamina
    (0, 0, 255),  # 3 = Green sheath
    (255, 255, 0),  # 4 = Senescent sheath
    (0, 255, 255),  # 5 = Green internode
    (255, 0, 255),  # 6 = senescent Internode
    (128, 255, 0),  # 7 = green peduncule
    (0, 128, 255),  # 8 = senescent peduncule
    (255, 0, 128),  # 9 = green ear
    (0, 255, 128),  # 10 = senescent ear
    (128, 0, 255),  # 11 = green awn
    (255, 128, 0),  # 12 = senescent awn
    (128, 128, 255),  # ???
    (255, 128, 128),
    (128, 255, 128),
    (255, 255, 255),
]


def col_item(ind, color_list=color_list):
    if ind is None:
        return lambda x: color_list[(x - 1) % len(color_list)]
    else:
        return (color_list[(ind - 1) % len(color_list)],)


def povray(
    scene,
    pov_file="./scene.pov",
    camera_distance=1.0,
    fov=45.0,
    width=320,
    height=280,
    domain=((-0.5, -0.5), (0.5, 0.5)),
    azimuth=0.0,
    zenith=0.0,
    camera_type="perspective",
    soil=False,
    povray_cmd="povray",
):
    """
     !!!! Deprecated function, use openalea.povray.Povray class instead !!!!

    Compute povray files based both on a scene and its stand box.

    :Parameters:
        - scene: a plantgl scene
        - camera distance: distance from the position of the camera to the look_at point
        - angle : angle corresponding to the width of the image
        - width: width of the final image in pixel
        - height: height of the final image in pixel
        - domain: scene pattern used by caribu
        - azimuth: angle in degree around the vertical axis. az=0. is equialent to have the width of the image align to X direction of the scene
        - zenith: angle between the view direction and the vertical
        - camera type: perspective, orthographic or fisheye
        - soil: add a soil to the scene
        - povray cmd : the path of the povray exe
    """

    f = path(pov_file)
    namebase = f.namebase
    ext = f.ext
    dirname = f.dirname()

    xc = 0.5 * (domain[0][0] + domain[1][0])
    yc = 0.5 * (domain[0][1] + domain[1][1])

    camera = {
        "type": camera_type,
        "distance": camera_distance,
        "fov": fov,
        "xc": xc,
        "yc": yc,
        "azimuth": azimuth,
        "zenith": zenith,
    }

    pov = PovRay(
        camera=camera, image_width=width, image_height=height, working_dir=str(dirname)
    )

    pov.render(scene, namebase + ext)
    image_name = pov.rendered_image_path

    d3D = domain3D(domain, scene)
    scene_box = pgl.Scene()
    scene_box.add(stand_box(d3D))
    f_box = namebase + "_box" + ext
    pov.render(scene_box, f_box)
    image_name_box = pov.rendered_image_path

    return image_name, image_name_box
