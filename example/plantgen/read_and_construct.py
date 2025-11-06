# -*- coding: latin-1 -*-
# -*- python -*-

"""
read_and_construct.py
~~~~~~~~~~~~~~~~~~~~~

An example to show how to construct the input data of Adel using :mod:`openalea.adel.plantgen`.

In this example, the inputs of PlantGen are read from the files 'inputs.py', 'dynT_user_MIN.csv'
and 'dimT_user_MIN.csv' using :mod:`openalea.adel.plantgen.read_plantgen_inputs`.

You must first install :mod:`openalea.adel` (and add it to your PYTHONPATH)
before running this script with the command `python read_and_construct.py`.

:copyright: Copyright 2015 INRA-EcoSys, Camille Chambon <camille.chambon@inra.fr>
:license: Cecill-C License, http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html

"""

# import the Pandas library. In this example, Pandas is used to read and
# write the tables.
import pandas as pd

# import the front end of PlanGen
from openalea.adel.plantgen import plantgen_interface

# read the inputs from 'inputs.py', 'dynT_user_MIN.csv' and 'dimT_user_MIN.csv'
(
    dynT_user,
    dimT_user,
    plants_number,
    plants_density,
    decide_child_axis_probabilities,
    MS_leaves_number_probabilities,
    ears_density,
    GL_number,
    delais_TT_stop_del_axis,
    TT_hs_break,
    inner_params,
) = plantgen_interface.read_plantgen_inputs(
    "inputs.py", "dynT_user_MIN.csv", "dimT_user_MIN.csv"
)

# run the construction of Adel inputs
axeT_, dimT_, phenT_, _, _, _, _, _, _, _ = plantgen_interface.gen_adel_input_data(
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
axeT_.to_csv("axeT.csv", na_rep="NA", index=False)
dimT_.to_csv("dimT.csv", na_rep="NA", index=False)
phenT_.to_csv("phenT.csv", na_rep="NA", index=False)

# "axeT.csv", "dimT.csv" and "phenT.csv" are now ready to be used by Adel.
