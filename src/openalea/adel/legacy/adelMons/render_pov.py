import os
from os.path import join

try:
    import win32api
except ImportError:
    win32api = None

from pathlib import Path as path
import time


class render_pov:
    """render a pov file and generates a bmp file . !! Pov-ray has to be launch !!"""

    def __call__(self, pov_file, povdir, wait):
        path_ini = os.path.abspath(os.getcwd())
        if win32api:
            dir = pov_file[
                : -len(pov_file.split("\\")[-1]) - 1
            ]  # path du fichier .pov d'origine
            name = pov_file.split("\\")[-1].split(".")[0]
            # os.system('copy '+win32api.GetShortPathName(pov_file)+' '+win32api.GetShortPathName(povdir))  #copie du fichier pov dans le bin de povray pas necessaire : avec short paths peut executer diretcement
            os.chdir(povdir)
            os.system(
                "pvengine.exe /RENDER "
                + join(win32api.GetShortPathName(dir), pov_file.split("\\")[-1])
            )  # short name pour le path, mais pas le nom de fichier pour conserver noms'd'image complets
            ext = ".bmp"
        else:
            _pov_file = path(pov_file)
            dir = _pov_file.dirname()
            name = _pov_file.namebase

            os.chdir(dir)
            ext = ".png"

        time.sleep(wait)  # temps d'attente en secondes entre 2 simulations
        os.chdir(path_ini)

        return join(dir, name + ext)
