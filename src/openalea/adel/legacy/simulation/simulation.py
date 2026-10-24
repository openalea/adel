import os
from pathlib import Path as path


def simulation(lsystem_filename, nb_plants):
    # write AleaPars.h
    fn = path(lsystem_filename)
    directory = fn.dirname()

    param_fn = directory / "AleaPars.h"
    s = "\n#define NB_PLANTS %d\n" % nb_plants
    f = open(param_fn, "w")
    f.write(s)
    f.close()

    view_file = "view.v"
    leaf_file = "leaf.a"
    color_file = "color.map"
    result_file = "result.str"
    cmd = "cpfg -g -homo -str %s %s %s %s -m %s" % (
        result_file,
        fn.name,
        view_file,
        leaf_file,
        color_file,
    )
    print(cmd)
    curdir = os.getcwd()
    print(directory)
    os.chdir(directory)
    os.system(cmd)
    os.chdir(curdir)
    result_file = directory / result_file
    if not result_file.exists():
        result_file = directory / "100plants" / "bid050.str"

    f = open(result_file)
    s = f.read()
    s = s.replace("\n", "")
    f.close()
    return s
