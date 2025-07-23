# -*- python -*-
#
#       Adel.PlantGen
#
#       Copyright 2012-2014 INRIA - CIRAD - INRA
#
#       File author(s): Camille Chambon <camille.chambon@grignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################
"""
This module stores the default values of :mod:`plantgen <openalea.adel.plantgen>`
inner parameters.

Authors: M. Abichou, B. Andrieu, C. Chambon
"""

SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS = {"a_1": 0.9423, "a_2": 0.555}
"""The coefficients *a_1* and *a_2* to calculate the final number of leaves on 
tillers from the final number of leaves on main stem. 

Calculation is done as follows::
    
    tiller_final_leaves_number 
        = a_1 * MS_final_leaves_number - a_2 * cohort_number
"""

MS_EMERGENCE_STANDARD_DEVIATION = 0.3
"""the standard deviation on the date of emergence of the main stems in the plot (in number of phyllochrons units).

This parameter is used to calculate the TT_*_phytomer1 value of the main stems."""

TILLERS_EMERGENCE_STANDARD_DEVIATION = 0.2
"""the standard deviation on the date of emergence of the tillers in the plot (in phyllochron units).

This parameter is used to calculate the TT_*_phytomer1 value of the tillers."""

MS_HS_AT_TILLER_EMERGENCE = {
    "T0": 1.9,
    "T1": 2.5,
    "T2": 3.2,
    "T3": 4.3,
    "T4": 5.5,
    "T5": 6.9,
    "T6": 9.16,
    "T7": 11.9,
    "T8": 15.4,
    "T9": 20.1,
    "T10": 26.1,
}
"""The haunstage of the most frequent main stem at the emergence of each tiller (in phyllochron units).

This parameter is a Python dictionary. 
The keys represent the cohort indexes and the values represent the phyllochronic delay. 
"""

N2_MS_DIV_N2_COHORT = 0.85
"""Ratio between the maximum number of green leaves on the tillers and the 
maximum green leaves on the main stem.

Value is given for the axes with the most frequent leaves number.
"""

DELAIS_PHYLL_COL_TIP_1ST = 1.0
"""Delay between tip appearance and collar appearance for the first leaf only.

The delay is given in phyllochron unit. 

To parameterize the delay between tip appearance and collar appearance for the other leaves, see 
:attr:`DELAIS_PHYLL_COL_TIP_NTH <openalea.adel.plantgen.params.DELAIS_PHYLL_COL_TIP_NTH>`.
"""

DELAIS_PHYLL_COL_TIP_NTH = 1.6
"""Delay between tip appearance and collar appearance for all leaves except the 
first one. 

The delay is given in phyllochron unit. 

To parameterize the delay between tip appearance and collar appearance for the first leaf, see 
:attr:`DELAIS_PHYLL_COL_TIP_1ST <openalea.adel.plantgen.params.DELAIS_PHYLL_COL_TIP_1ST>`.   
"""

DELAIS_PHYLL_HS_COL_NTH = -0.20
"""Delay between Haun Stage and collar appearance for all leaves. 

The delay is given in phyllochron unit. 
"""

DELAIS_PHYLL_SEN_DISP = 3.0
"""The time during which a fully senesced leaf on a non-elongated internode 
remains on the plant. 

The delay is given in phyllochron unit. 
"""

TT_DEL_FHAUT = 3000
"""The thermal time at which leaves on elongated internode disappear.

The thermal time is given in degree.day. 
"""


FIRST_CHILD_DELAY = 2
"""The delay between a parent cohort and its first possible child cohort. 
This delay is expressed in number of cohorts."""

LENGTHS_REDUCTION_FACTOR = -0.115
"""The reduction factor to apply on the lengths of the most frequent main stem to compute the lengths of the regressive tillers,
i.e. lengths_regressive_tillers =  lengths_most_frequent_MS * (1 + LENGTHS_REDUCTION_FACTOR)
"""

WIDTHS_REDUCTION_FACTOR = -0.0938
"""The reduction factor to apply on the widths of the most frequent main stem to compute the widths of the regressive tillers,
i.e. widths_regressive_tillers =  widths_most_frequent_MS * (1 + WIDTHS_REDUCTION_FACTOR)
"""

SLOPE_SHIFT_MS_TO_TILLERS = -0.6681
"""Shift to apply on phytomer index for each cohort.
"""

TILLERS_L_BLADE_1ST = 6.5
"""The length of the blades of phytomer 1 for every tiller.
"""

W_INTERNODE_POLYNOMIAL = {
    "coefficients": {"a0": -2.4527, "a1": 7.6398, "a2": -4.207},
    "first_point": {
        "index_relative_to_MS_phytomer_normalized": 0.6,
        "W_internode_normalized": 0.6,
    },
}
"""Parameters used to compute the width of the internodes for each tiller:
    - the coefficients of the polynomial fitted from main stems experimental internode widths,
    - the first point where the polynomial can be applied."""

MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY = 0.0  # 0.26
"""The delay to apply on SSI of the most frequent MS to compute SSI of regressive tillers, 
i.e. SSI_regressive_tillers = SSI_most_frequent_MS + MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
"""

A_COHORT_REDUCTION_FACTOR = -0.53
"""The reduction factor to apply on the a_cohort of the most frequent main stem to compute the a_cohort of the regressive tillers,
i.e. a_cohort_regressive_tillers =  a_cohort_most_frequent_MS * (1 + A_COHORT_REDUCTION_FACTOR)
"""

FLAG_LIGULATION_DELAY = 20
"""The ligulation delay of the flag. Used to compute dTT_MS_cohort.
"""

EMERGENCE_PROBABILITY_REDUCTION_FACTOR = 0.29
"""The reduction factor of the emergence probability of secondary tiller compared to primary one,
i.e. probability(secondary tiller)= (1 - EMERGENCE_PROBABILITY_REDUCTION_FACTOR) * probability(primary tiller).
e.g. : probability(T1.0.0) = probability(T1) * (1 - EMERGENCE_PROBABILITY_REDUCTION_FACTOR) * probability(T3) * (1 - EMERGENCE_PROBABILITY_REDUCTION_FACTOR) * probability(T5) 
"""

NUMBER_OF_ELONGATED_INTERNODES = 4
"""The number of elongated internodes. 
This number is used to compute the time (t1) corresponding to the minimum of green leaves number(n1) for the main stem and tillers.
"""

K1 = 0.854
"""Proportional factor of blade width relative to the first leaf of each tiller.
"""

K2 = 0.91
"""Proportional factor of blade width relative to the last leaf of each tiller.
"""

START_MS_HS_MORTALITY_VS_N_PHYTOMER = {"as": 0.2411, "bs": 4.2546}
"""Correlation between the start of Haun stage mortality of the main stem and N_phytomer of the main stem.
"""

END_MS_HS_MORTALITY_VS_N_PHYTOMER = {"ae": 1.0247}
"""Correlation between the end of Haun stage mortality of the main stem and N_phytomer of the main stem.
"""
