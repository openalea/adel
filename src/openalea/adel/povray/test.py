from openalea.adel.povray.povray import PovRay, povray
import openalea.plantgl.all as pgl
import cv2


# pov = PovRay(working_dir='./test')


def test_class_interface():
    scene = pgl.Scene()
    scene.add(pgl.Sphere(radius=0.1))
    pov = PovRay()
    pov.render(scene)
    reader = cv2.imread
    im = pov.get_image(reader)
    return im


def test_old_interface():
    scene = pgl.Scene()
    scene.add(pgl.Sphere(radius=0.1))
    return povray(scene)
