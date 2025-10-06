# -*- coding: latin-1 -*-
# -*- python -*-

"""
construct.py
~~~~~~~~~~~~

An example to show how to construct the input data of Adel using :mod:`openalea.adel.plantgen`.

In this example, the inputs of PlantGen are hard coded in the code of the example.
See the script read_and_construct.py for an example of how to read the inputs of PlantGen from a file.

You must first install :mod:`openalea.adel` (and add it to your PYTHONPATH)
before running this script with the command `python construct.py`.

:copyright: Copyright 2015 INRA-EcoSys, Camille Chambon <camille.chambon@grignon.inra.fr>
:license: Cecill-C License, http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html

"""

# import the Pandas library. In this example, Pandas is used to read and
# write the tables.
import pandas as pd

# read the dynT_user_MIN table. "dynT_user_MIN.csv" must be in the working directory.
dynT_user = pd.read_csv("dynT_user_MIN.csv")

# read the dimT_user_MIN table. "dimT_user_MIN.csv" must be in the working directory.
dimT_user = pd.read_csv("dimT_user_MIN.csv")

# define the other inputs
plants_number = 100
plants_density = 250
decide_child_axis_probabilities = {
    "T0": 0.0,
    "T1": 0.900,
    "T2": 0.983,
    "T3": 0.817,
    "T4": 0.117,
}
MS_leaves_number_probabilities = {
    "10": 0.145,
    "11": 0.818,
    "12": 0.037,
    "13": 0.0,
    "14": 0.0,
}
ears_density = 500
GL_number = {1117.0: 5.6, 1212.1: 5.4, 1368.7: 4.9, 1686.8: 2.4, 1880.0: 0.0}
delais_TT_stop_del_axis = 600
inner_params = {"DELAIS_PHYLL_COL_TIP_1ST": 1.0, "DELAIS_PHYLL_COL_TIP_NTH": 1.6}

# run the construction of Adel inputs
from openalea.adel.plantgen import plantgen_interface

axeT, dimT, phenT, _, _, _, _, _, _, _ = plantgen_interface.gen_adel_input_data(
    dynT_user,
    dimT_user,
    plants_number,
    plants_density,
    decide_child_axis_probabilities,
    MS_leaves_number_probabilities,
    ears_density,
    GL_number,
    delais_TT_stop_del_axis,
    inner_params=inner_params,
)

# write axeT, dimT and phenT to CSV files in the working directory, replacing
# missing values by 'NA' and ignoring the indexes
axeT.to_csv("axeT.csv", na_rep="NA", index=False)
dimT.to_csv("dimT.csv", na_rep="NA", index=False)
phenT.to_csv("phenT.csv", na_rep="NA", index=False)

# "axeT.csv", "dimT.csv" and "phenT.csv" are now ready to be used by Adel.
