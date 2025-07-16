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
Routines defining the main steps of the process.

Routines with a leading underscore are non-public and shouldn't be used by the
user. They defines low-level processes, and must be called in a specific order.

Some of the routines take an optional "force" input parameter. By default this parameter
has the value True.
When "force" is True, then the variables returned by the routine are recalculated.
When "force" is False, the routine first checks if the variables have already been g
calculated before. If so, the variables are not recalculated.
The parameter "force" primary aims to not recalculate the random variables, so to
avoid incoherent results from different calls.
The parameter "force" also permits to avoid recalculation when it is useless.

Authors: M. Abichou, B. Andrieu, C. Chambon
"""

import math
import random

import numpy as np
import pandas as pd

from alinea.adel.plantgen import tools, params


class DataCompleteness:
    """
    This enumerate defines the different degrees of completeness that the data
    documented by the user can have.
    """

    MIN = "MIN"
    SHORT = "SHORT"
    FULL = "FULL"


def init_axes(
    plants_number,
    decide_child_cohort_probabilities,
    MS_leaves_number_probabilities,
    theoretical_cohort_cardinalities,
    theoretical_axis_cardinalities,
    axeT_user=None,
):
    """
    Initialize the axes randomly.
    The following variables are calculated:

        * id_plt: Number (int) identifying the plant to which the axe belongs
        * id_cohort: Number (int) identifying the cohort to which the axe belongs
        * id_axis: Identifier of the botanical position of the axis on the plant.
              MS refers to the main stem, T0, T1, T2,..., refers to the primary tillers,
              T0.0, T0.1, T0.2,..., refers to the secondary tillers of the primary tiller
              T0, and T0.0.0, T0.0.1, T0.0.2,..., refers to the tertiary tillers
              of the secondary tiller T0.0.
        * N_phytomer_potential: The potential total number of vegetative phytomers formed on
            the axis. N_phytomer_potential does NOT take account of the regression of some axes.
        * id_phen: a key (int) linking to phenT. id_phen allows referring to the data
            that describe the phenology of the axis

    and are stored in memory for the next steps of the process.
    The routine returns :ref:`cardinalityT` for debugging purpose.
    """
    # create axeT_tmp
    if axeT_user is None:
        axeT_tmp = _create_axeT_tmp(
            plants_number,
            decide_child_cohort_probabilities,
            MS_leaves_number_probabilities,
        )
    else:
        axeT_tmp = axeT_user
        _create_axeT_tmp.axeT_tmp = axeT_user

    # create cardinalityT
    cardinalityT = _create_cardinalityT(
        theoretical_cohort_cardinalities,
        theoretical_axis_cardinalities,
        axeT_tmp[["id_cohort", "id_axis"]],
    )

    return cardinalityT


class PhenologyFunctions:
    """
    Define phenology functions. These functions are used later on to calculate the
    phenology of the axes.
    The following variables are calculated:

        * cardinality: the cardinality of the couple (id_axis, N_phytomer_potential) in axeT
        * a_cohort: the rate of Haun Stage vs Thermal time. This is the rate of the first phase in case of bilinear behavior.
        * TT_hs_0: the thermal time for Haun Stage equal to 0
        * TT_hs_break: the thermal time when the rate of phytomers emergence is changing. NA: no change.
        * TT_flag_ligulation: the thermal time when Haun Stage is equal to N_phytomer_potential
        * dTT_MS_cohort: the delays between the emergence of the main stem and the emergence of the cohort.
        * n0: number of green leaves at t0
        * n1: number of green leaves at t1
        * n2: number of green leaves at TT_flag_ligulation
        * t0: the thermal time at the start of leaf senescence
        * t1: the thermal time at which the senescence starts
        * hs_t1: the Haun Stage at t1
        * a: the coefficient of the 3rd order term of the polynomial describing the dynamics
            of the number of green leaves after flowering
        * c: the coefficient of the 1st order term of the polynomial describing the dynamics
            of the number of green leaves after flowering
        * RMSE_gl: the RMSE for the dynamic of the number of green leaves after estimation of
            parameter a.

    and are stored in memory for the next steps of the process.
    The routine returns :ref:`dynT` and the decimal number of elongated internodes.
    """

    def __init__(self):
        self.dynT_ = None
        self.decimal_elongated_internode_number = None

    def __call__(
        self,
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        dynT_user_completeness,
        dimT_user_completeness,
        TT_hs_break,
        force=True,
        axeT_user=None,
        TT_t1_user=None,
    ):
        if force or self.dynT_ is None:
            # 1. create axeT_tmp, dynT_tmp and dimT_tmp
            if axeT_user is None:
                axeT_tmp = _create_axeT_tmp(
                    plants_number,
                    decide_child_cohort_probabilities,
                    MS_leaves_number_probabilities,
                    force=False,
                )
            else:
                axeT_tmp = axeT_user
            dynT_tmp = _create_dynT_tmp(axeT_tmp)

            dimT_tmp = _create_dimT_tmp(axeT_tmp)

            # 2. merge dynT_tmp and dynT_user
            dynT_tmp_merged = _merge_dynT_tmp_and_dynT_user(
                dynT_tmp, dynT_user, dynT_user_completeness, TT_hs_break
            )

            # 3. merge dimT_tmp and dimT_user
            dimT_tmp_merged = _merge_dimT_tmp_and_dimT_user(
                dynT_tmp_merged, dimT_user, dimT_user_completeness, dimT_tmp
            )

            # 4. calculate decimal_elongated_internode_number
            self.decimal_elongated_internode_number = (
                _calculate_decimal_elongated_internode_number(
                    dimT_tmp_merged, dynT_tmp_merged
                )
            )

            # 5. create dynT
            self.dynT_ = _create_dynT(dynT_tmp_merged, GL_number, TT_t1_user=TT_t1_user)

        return self.dynT_, self.decimal_elongated_internode_number


phenology_functions = PhenologyFunctions()


def plants_structure(
    plants_number,
    decide_child_cohort_probabilities,
    MS_leaves_number_probabilities,
    dynT_user,
    dimT_user,
    GL_number,
    dynT_user_completeness,
    dimT_user_completeness,
    TT_hs_break,
    delais_TT_stop_del_axis,
    number_of_ears,
    plants_density,
    ears_density,
    axeT_user=None,
    TT_regression_start_user=None,
    TT_t1_user=None,
):
    """
    Construct the structure of the plants.
    The following variables are calculated:

        * N_phytomer: The effective total number of vegetative phytomers formed
            on the axis. N_phytomer does take account of the regression of some axes.
        * HS_final: The Haun Stage at the end of growth of the axis.
        * TT_stop_axis: If the axis dyes: thermal time (since crop emergence) of
            end of growth. If the axis grows up to flowering: NA
        * TT_del_axis: If the axis dyes: thermal time (since crop emergence) of
            disappearance. If the axis grows up to flowering: NA
        * id_dim: key (int) linking to dimT. id_dim allows referring to the data
            that describe the dimensions of the phytomers of the axis
        * id_ear: Key (int) linking to earT. id_ear allows referring to the data
            that describe the ear of the axis. For the regressive axes, id_ear=NA.
            For the non-regressive axes, id_ear=1.
        * TT_em_phytomer1: Thermal time (relative to canopy appearance) of tip
            appearance of the first true leaf (not coleoptile or prophyll)
        * TT_col_phytomer1: Thermal time (relative to canopy appearance) of collar
            appearance of the first true leaf
        * TT_sen_phytomer1: Thermal time (relative to canopy appearance) of full
            senescence of the first true leaf (this is : thermal time when SSI= 1)
        * TT_del_phytomer1: Thermal time (relative to canopy appearance) of
            disappearance of the first true leaf

    and are stored in memory for the next steps of the process.
    The routine returns

        * :ref:`axeT <axeT>` as final result,
        * :ref:`tilleringT` and :ref:`phenT_first` for debugging purpose.
    """
    # 1. create axeT_tmp, dynT_
    if axeT_user is None:
        axeT_tmp = _create_axeT_tmp(
            plants_number,
            decide_child_cohort_probabilities,
            MS_leaves_number_probabilities,
            force=False,
        )
    else:
        axeT_tmp = axeT_user

    dynT_, decimal_elongated_internode_number = phenology_functions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        dynT_user_completeness,
        dimT_user_completeness,
        TT_hs_break,
        force=False,
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )

    # 2. create phenT_tmp
    phenT_tmp = _create_phenT_tmp(axeT_tmp, dynT_, decimal_elongated_internode_number)

    # 3. create phenT_first
    phenT_first = _create_phenT_first(phenT_tmp)

    # 4. create axeT
    axeT_, axeT_tmp, phenT_tmp, phenT_first, TT_regression_start, TT_regression_end = (
        _create_axeT(
            axeT_tmp,
            phenT_first,
            dynT_,
            delais_TT_stop_del_axis,
            number_of_ears,
            TT_regression_start_user=TT_regression_start_user,
        )
    )

    # 5. create tilleringT
    tilleringT = _create_tilleringT(
        dynT_,
        phenT_first,
        axeT_.index.size,
        plants_number,
        plants_density,
        ears_density,
        TT_regression_start,
        TT_regression_end,
    )

    return axeT_, tilleringT, phenT_first


def organs_dimensions(
    plants_number,
    decide_child_cohort_probabilities,
    MS_leaves_number_probabilities,
    dynT_user,
    dimT_user,
    GL_number,
    dynT_user_completeness,
    dimT_user_completeness,
    TT_hs_break,
    delais_TT_stop_del_axis,
    number_of_ears,
    axeT_user=None,
    TT_t1_user=None,
):
    """
    Calculate the dimensions of the organs.
    The following variables are calculated:

        * id_dim: key (int) of the axis
        * index_phytomer: The absolute phytomer position, i.e. phytomer rank
        * L_blade: length of the mature blade (cm)
        * W_blade: Maximum width of the mature leaf blade (cm)
        * L_sheath: Length of a mature sheath (cm)
        * W_sheath: Diameter of the stem or pseudo stem at the level of sheath (cm)
        * L_internode: Length of an internode (cm)
        * W_internode: Diameter of an internode (cm)

    and are stored in memory for the next steps of the process.
    The routine returns :ref:`dimT <dimT>` as final result.
    """
    # 1. create axeT_, dimT_tmp, phenT_tmp and dynT_
    if axeT_user is None:
        axeT_tmp = _create_axeT_tmp(
            plants_number,
            decide_child_cohort_probabilities,
            MS_leaves_number_probabilities,
            force=False,
        )
    else:
        axeT_tmp = axeT_user
    dynT_, decimal_elongated_internode_number = phenology_functions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        dynT_user_completeness,
        dimT_user_completeness,
        TT_hs_break,
        force=False,
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )
    phenT_tmp = _create_phenT_tmp(
        axeT_tmp, dynT_, decimal_elongated_internode_number, force=False
    )
    phenT_first = _create_phenT_first(phenT_tmp, force=False)
    axeT_, axeT_tmp, phenT_tmp, phenT_first, TT_regression_start, TT_regression_end = (
        _create_axeT(
            axeT_tmp,
            phenT_first,
            dynT_,
            delais_TT_stop_del_axis,
            number_of_ears,
            force=False,
        )
    )
    dimT_tmp = _create_dimT_tmp(axeT_tmp, force=False)

    # 2. create dimT
    dimT_ = _create_dimT(axeT_, dimT_tmp, dynT_, decimal_elongated_internode_number)

    return dimT_


def axes_phenology(
    plants_number,
    decide_child_cohort_probabilities,
    MS_leaves_number_probabilities,
    dynT_user,
    dimT_user,
    GL_number,
    dynT_user_completeness,
    dimT_user_completeness,
    TT_hs_break,
    delais_TT_stop_del_axis,
    number_of_ears,
    axeT_user=None,
    TT_t1_user=None,
):
    """
    Calculate the phenology of the axes.
    The following variables are calculated:

        * TT_em_phytomer: Thermal time of the appearance of the tip of leaf out of
            the whorl made by the older blade.
        * TT_col_phytomer: Thermal time of the appearance of collar.
        * TT_sen_phytomer: Thermal time for which SSI = n (where n is the phytomer
             rank).
        * TT_del_phytomer: Thermal time after which the leaf blade is destroyed
                and is not displayed in the 3D mock-up anymore.
        * dTT_em_phytomer: Thermal time of the appearance of the tip of leaf out of
                the whorl made by the older blade; expressed as thermal time since TT_em_phytomer1
        * dTT_col_phytomer: Thermal time of the appearance of collar; expressed as
                thermal time since TT_col_phytomer1
        * dTT_sen_phytomer: Thermal time for which SSI = n (where n is the phytomer
                rank); expressed as thermal time since TT_sen_phytomer1
        * dTT_del_phytomer: Thermal time after which the leaf blade is destroyed
              and is not displayed in the 3D mock-up anymore; expressed as thermal time
              since TT_del_phytomer1

    and are stored in memory for the next steps of the process.
    The routine returns:

        * :ref:`phenT <phenT_>` as final result,
        * :ref:`phenT_abs` and :ref:`HS_GL_SSI_T` for debugging purpose.
    """
    # 1. create phenT_tmp, axeT_, dimT_, phenT_first and dynT_
    if axeT_user is None:
        axeT_tmp = _create_axeT_tmp(
            plants_number,
            decide_child_cohort_probabilities,
            MS_leaves_number_probabilities,
            force=False,
        )
    else:
        axeT_tmp = axeT_user
    dynT_, decimal_elongated_internode_number = phenology_functions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        dynT_user_completeness,
        dimT_user_completeness,
        TT_hs_break,
        force=False,
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )
    phenT_tmp = _create_phenT_tmp(
        axeT_tmp, dynT_, decimal_elongated_internode_number, force=False
    )
    phenT_first = _create_phenT_first(phenT_tmp, force=False)
    axeT_, axeT_tmp, phenT_tmp, phenT_first, TT_regression_start, TT_regression_end = (
        _create_axeT(
            axeT_tmp,
            phenT_first,
            dynT_,
            delais_TT_stop_del_axis,
            number_of_ears,
            force=False,
        )
    )
    dimT_tmp = _create_dimT_tmp(axeT_tmp, force=False)
    dimT_ = _create_dimT(
        axeT_, dimT_tmp, dynT_, decimal_elongated_internode_number, force=False
    )

    # create phenT_abs
    phenT_abs = _create_phenT_abs(phenT_tmp, axeT_, dimT_)
    # create phenT
    phenT_ = _create_phenT(phenT_abs, phenT_first)
    # create HS_GL_SSI_T
    HS_GL_SSI_T = _create_HS_GL_SSI_T(axeT_, dynT_)

    return phenT_, phenT_abs, HS_GL_SSI_T


################################################################################

############### PRIVATE : DO NOT USE FROM OUTSIDE THE MODULE ###################

################################################################################


class _CreateAxeTTmp:
    """
    Create the *axeT_tmp* dataframe.
    Compute the following columns: *id_plt*, *id_cohort*, *id_axis*, *N_phytomer_potential* and *id_phen*.
    """

    def __init__(self):
        self.axeT_tmp = None

    def __call__(
        self,
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        force=True,
    ):
        if force or self.axeT_tmp is None:
            plant_ids = list(range(1, plants_number + 1))
            id_cohort_list, id_axis_list = _gen_id_axis_list(
                plant_ids, decide_child_cohort_probabilities
            )
            id_plt_list = _gen_id_plt_list(plant_ids, id_cohort_list)
            N_phytomer_potential_list = _gen_N_phytomer_potential_list(
                id_cohort_list,
                MS_leaves_number_probabilities,
                params.SECONDARY_STEM_LEAVES_NUMBER_COEFFICIENTS,
            )
            id_phen_list = _gen_id_phen_list(id_cohort_list, N_phytomer_potential_list)

            self.axeT_tmp = pd.DataFrame(
                index=list(range(len(id_plt_list))),
                columns=[
                    "id_plt",
                    "id_cohort",
                    "id_axis",
                    "N_phytomer_potential",
                    "N_phytomer",
                    "HS_final",
                    "TT_stop_axis",
                    "TT_del_axis",
                    "id_dim",
                    "id_phen",
                    "id_ear",
                    "TT_em_phytomer1",
                    "TT_col_phytomer1",
                    "TT_sen_phytomer1",
                    "TT_del_phytomer1",
                ],
                dtype=float,
            )
            self.axeT_tmp = self.axeT_tmp.astype({'id_axis':'str'})
            self.axeT_tmp.loc[:,"id_plt"] = id_plt_list
            self.axeT_tmp.loc[:,"id_cohort"] = id_cohort_list
            self.axeT_tmp.loc[:,"id_axis"] = id_axis_list
            self.axeT_tmp.loc[:,"N_phytomer_potential"] = N_phytomer_potential_list
            self.axeT_tmp.loc[:,"N_phytomer"] = N_phytomer_potential_list
            self.axeT_tmp.loc[:,"id_phen"] = id_phen_list
        return self.axeT_tmp


_create_axeT_tmp = _CreateAxeTTmp()


class _CreateAxeT:
    """
    Create the :ref:`axeT <axeT>` dataframe, and update axeT_tmp, phenT_tmp and phenT_first according to regressive tillers.
    """

    def __init__(self):
        self.axeT_ = None
        self.TT_regression_start = None
        self.TT_regression_end = None

    def __call__(
        self,
        axeT_tmp,
        phenT_first,
        dynT_,
        delais_TT_stop_del_axis,
        number_of_ears,
        force=True,
        TT_regression_start_user=None,
    ):
        if (
            force
            or self.axeT_ is None
            or self.TT_regression_start is None
            or self.TT_regression_end is None
        ):
            self.axeT_ = axeT_tmp.copy()

            TT_hs_break, N_phytomer_potential, a_cohort, TT_hs_0, TT_flag_ligulation = (
                dynT_.loc[
                    dynT_.first_valid_index(),
                    [
                        "TT_hs_break",
                        "N_phytomer_potential",
                        "a_cohort",
                        "TT_hs_0",
                        "TT_flag_ligulation",
                    ],
                ]
            )

            if TT_regression_start_user is None:
                as_ = params.START_MS_HS_MORTALITY_VS_N_PHYTOMER["as"]
                bs_ = params.START_MS_HS_MORTALITY_VS_N_PHYTOMER["bs"]
                MS_HS_tillering_mortality_start = as_ * N_phytomer_potential + bs_
                if math.isnan(TT_hs_break):  # linear mode
                    self.TT_regression_start = (
                        MS_HS_tillering_mortality_start / a_cohort + TT_hs_0
                    )
                else:  # bilinear mode
                    HS_break_0 = a_cohort * (TT_hs_break - TT_hs_0)
                    if MS_HS_tillering_mortality_start < HS_break_0:  # 1rst phase
                        self.TT_regression_start = (
                            MS_HS_tillering_mortality_start / a_cohort + TT_hs_0
                        )
                    else:  # 2nd phase
                        a2_0 = (N_phytomer_potential - HS_break_0) / (
                            TT_flag_ligulation - TT_hs_break
                        )
                        self.TT_regression_start = (
                            MS_HS_tillering_mortality_start - HS_break_0
                        ) / a2_0 + TT_hs_break
            else:
                self.TT_regression_start = TT_regression_start_user

            ae_ = params.END_MS_HS_MORTALITY_VS_N_PHYTOMER["ae"]
            MS_HS_tillering_mortality_end = ae_ * N_phytomer_potential
            if math.isnan(TT_hs_break):  # linear mode
                self.TT_regression_end = (
                    MS_HS_tillering_mortality_end / a_cohort + TT_hs_0
                )
            else:  # bilinear mode
                HS_break_0 = a_cohort * (TT_hs_break - TT_hs_0)
                if MS_HS_tillering_mortality_end < HS_break_0:  # 1rst phase
                    self.TT_regression_end = (
                        MS_HS_tillering_mortality_end / a_cohort + TT_hs_0
                    )
                else:  # 2nd phase
                    a2_0 = (N_phytomer_potential - HS_break_0) / (
                        TT_flag_ligulation - TT_hs_break
                    )
                    self.TT_regression_end = (
                        MS_HS_tillering_mortality_end - HS_break_0
                    ) / a2_0 + TT_hs_break

            (
                self.axeT_.loc[:,"TT_em_phytomer1"],
                self.axeT_.loc[:,"TT_col_phytomer1"],
                self.axeT_.loc[:,"TT_sen_phytomer1"],
                self.axeT_.loc[:,"TT_del_phytomer1"],
            ) = _gen_all_TT_phytomer1_list(axeT_tmp, phenT_first, dynT_)
            self.axeT_.loc[self.axeT_["id_axis"] == "MS", "TT_stop_axis"] = np.nan
            tillers_axeT_index = self.axeT_.loc[self.axeT_["id_axis"] != "MS"].index
            self.axeT_.loc[tillers_axeT_index, "TT_stop_axis"] = (
                tools.decide_time_of_death(
                    axeT_tmp.index.size,
                    number_of_ears,
                    self.TT_regression_start,
                    self.TT_regression_end,
                    self.axeT_.loc[
                        tillers_axeT_index, ["id_plt", "id_axis", "TT_em_phytomer1"]
                    ],
                )
            )
            self.axeT_.loc[:,"id_ear"] = _gen_id_ear_list(self.axeT_["TT_stop_axis"])
            self.axeT_.loc[:,"TT_del_axis"] = _gen_TT_del_axis_list(
                self.axeT_["TT_stop_axis"], delais_TT_stop_del_axis
            )
            HS_final_series = _gen_HS_final_series(self.axeT_, dynT_)
            self.axeT_ = _remove_axes_without_leaf(self.axeT_, HS_final_series.index)
            self.axeT_.loc[:,"HS_final"] = HS_final_series.values
            self.axeT_.loc[:,"N_phytomer"] = _gen_N_phytomer(self.axeT_["HS_final"])
            self.axeT_.loc[:,"id_dim"] = _gen_id_dim_list(
                self.axeT_["id_cohort"], self.axeT_["N_phytomer"], self.axeT_["id_ear"]
            )
            self.axeT_.loc[:,"id_phen"] = self.axeT_["id_dim"]

            # we now know which tillers are regressive and which are not ; update axeT_tmp, phenT_tmp and phenT_first accordingly.
            _create_axeT_tmp.axeT_tmp.N_phytomer = self.axeT_.N_phytomer
            _create_axeT_tmp.axeT_tmp.id_phen = self.axeT_.id_phen
            phenT_tmp = _create_phenT_tmp(
                self.axeT_,
                dynT_,
                phenology_functions.decimal_elongated_internode_number,
            )
            _create_phenT_first(phenT_tmp)

        return (
            self.axeT_,
            _create_axeT_tmp.axeT_tmp,
            _create_phenT_tmp.phenT_tmp,
            _create_phenT_first.phenT_first,
            self.TT_regression_start,
            self.TT_regression_end,
        )


_create_axeT = _CreateAxeT()


def _gen_id_plt_list(plant_ids, id_cohort_list):
    """Generate the *id_plt* column."""
    id_plt_list = []
    current_plant_index = 0
    for plant_id in plant_ids:
        start_index = current_plant_index + 1
        if 1 in id_cohort_list[start_index:]:
            next_plant_first_row = id_cohort_list.index(1, start_index)
        else:
            next_plant_first_row = len(id_cohort_list)
        current_plant_axes = id_cohort_list[current_plant_index:next_plant_first_row]
        id_plt_list.extend([plant_id for _ in current_plant_axes])
        current_plant_index = next_plant_first_row
    return id_plt_list


def _gen_id_axis_list(plant_ids, decide_child_cohort_probabilities):
    """Generate the columns *id_axis* and *id_cohort* ."""
    all_child_cohorts = []
    for plant_id in plant_ids:
        child_cohorts = tools.decide_child_cohorts(
            decide_child_cohort_probabilities,
            params.FIRST_CHILD_DELAY,
            params.EMERGENCE_PROBABILITY_REDUCTION_FACTOR,
        )
        child_cohorts.sort()
        all_child_cohorts.extend(child_cohorts)
    all_child_cohorts_array = np.array(all_child_cohorts)
    cohort_numbers = all_child_cohorts_array[:, 0].astype(int).tolist()
    cohort_positions = all_child_cohorts_array[:, 1].tolist()
    return cohort_numbers, cohort_positions


def _gen_N_phytomer_potential_list(
    id_cohort_list,
    MS_leaves_number_probabilities,
    secondary_stem_leaves_number_coefficients,
):
    """Generate the *N_phytomer_potential* column."""
    N_phytomer_potential_list = []
    MS_final_leaves_number = 0.0
    # for each plant...
    for cohort_number in id_cohort_list:
        # calculate the leaves number of each axis
        leaves_number_float = 0.0
        if cohort_number == 1:
            # It is the main stem, then the leaves number has to satisfy the probability distribution defined
            # in MS_leaves_number_probabilities
            MS_final_leaves_number = tools.calculate_MS_final_leaves_number(
                MS_leaves_number_probabilities
            )
            leaves_number_float = MS_final_leaves_number
        else:
            # it is a secondary stem (i.e. a tiller)
            leaves_number_float = tools.calculate_tiller_final_leaves_number(
                MS_final_leaves_number,
                cohort_number,
                secondary_stem_leaves_number_coefficients,
            )
        fractional_part, integer_part = math.modf(leaves_number_float)
        if random.random() <= fractional_part:
            leaves_number_int = int(math.ceil(leaves_number_float))
        else:
            leaves_number_int = int(integer_part)
        N_phytomer_potential_list.append(leaves_number_int)

    return N_phytomer_potential_list


def _gen_N_phytomer(HS_final_series):
    """Generate the *N_phytomer* column."""
    return np.ceil(HS_final_series).astype(int)


def _gen_all_TT_phytomer1_list(axeT_tmp, phenT_first, dynT):
    """Generate the *TT_em_phytomer1*, *TT_col_phytomer1*, *TT_sen_phytomer1* and *TT_del_phytomer1* columns.
    For each plant, define a delay of appearance, and for each axis add this delay to the first leaf development schedule."""
    MS_sigma = params.MS_EMERGENCE_STANDARD_DEVIATION
    MS_sigma_div_2 = MS_sigma / 2.0
    tillers_sigma = params.TILLERS_EMERGENCE_STANDARD_DEVIATION
    tillers_sigma_div_2 = tillers_sigma / 2.0

    TT_em_phytomer1_series = pd.Series(index=axeT_tmp.index)
    TT_col_phytomer1_series = pd.Series(index=axeT_tmp.index)
    TT_sen_phytomer1_series = pd.Series(index=axeT_tmp.index)
    TT_del_phytomer1_series = pd.Series(index=axeT_tmp.index)

    for id_plt, axeT_tmp_grouped_by_id_plt in axeT_tmp.groupby("id_plt"):
        for (
            id_phen,
            N_phytomer_potential,
            id_cohort,
        ), axeT_tmp_grouped_by_id_plt_and_id_phen in axeT_tmp_grouped_by_id_plt.groupby(
            ["id_phen", "N_phytomer_potential", "id_cohort"]
        ):
            primary_axis = tools.get_primary_axis(
                max(axeT_tmp_grouped_by_id_plt_and_id_phen.id_axis),
                params.FIRST_CHILD_DELAY,
            )

            if primary_axis == "MS":
                sigma = MS_sigma
                sigma_div_2 = MS_sigma_div_2
            else:  # tillers
                sigma = tillers_sigma
                sigma_div_2 = tillers_sigma_div_2

            infimum = -sigma_div_2
            supremum = sigma_div_2

            normal_distribution = random.normalvariate(mu=0.0, sigma=sigma)
            while normal_distribution < infimum or normal_distribution > supremum:
                normal_distribution = random.normalvariate(mu=0.0, sigma=sigma)

            dynT_group = dynT.loc[
                (dynT.id_cohort == id_cohort)
                & (dynT.N_phytomer_potential == N_phytomer_potential)
            ]
            a_cohort = dynT_group.loc[dynT_group.first_valid_index(), "a_cohort"]
            normal_distribution_in_growing_degree_days = normal_distribution / a_cohort

            current_row = phenT_first[phenT_first["id_phen"] == id_phen]
            first_valid_index = current_row.first_valid_index()
            TT_em_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = (
                normal_distribution_in_growing_degree_days
                + current_row["TT_em_phytomer"][first_valid_index]
            )
            TT_col_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = (
                normal_distribution_in_growing_degree_days
                + current_row["TT_col_phytomer"][first_valid_index]
            )
            TT_sen_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = (
                normal_distribution_in_growing_degree_days
                + current_row["TT_sen_phytomer"][first_valid_index]
            )
            TT_del_phytomer1_series[axeT_tmp_grouped_by_id_plt_and_id_phen.index] = (
                normal_distribution_in_growing_degree_days
                + current_row["TT_del_phytomer"][first_valid_index]
            )

    return (
        TT_em_phytomer1_series,
        TT_col_phytomer1_series,
        TT_sen_phytomer1_series,
        TT_del_phytomer1_series,
    )


def _gen_id_dim_list(id_cohort_series, N_phytomer_series, id_ear_series):
    """Generate the *id_dim* column."""
    is_ear = pd.Series(0, index=id_ear_series.index)
    is_ear[id_ear_series.dropna().index] = 1
    zfilled_array = np.core.defchararray.zfill(np.char.mod("%d", N_phytomer_series), 2)
    id_cohort_str_array = np.char.mod("%d", id_cohort_series)
    id_dim_array = np.core.defchararray.add(id_cohort_str_array, zfilled_array)
    id_dim_array = np.core.defchararray.add(
        id_dim_array, np.char.mod("%d", is_ear)
    ).astype(int)
    return id_dim_array.tolist()


def _gen_id_phen_list(id_cohort_list, N_phytomer_potential_list):
    """Generate the *id_phen* column."""
    id_phen_list = []
    for i in range(len(id_cohort_list)):
        id_phen_list.append(
            int(
                "".join(
                    [
                        str(id_cohort_list[i]),
                        str(N_phytomer_potential_list[i]).zfill(2),
                        "1",
                    ]
                )
            )
        )  # 1: axis with ear
    return id_phen_list


def _gen_id_ear_list(TT_stop_axis):
    """Generate the *id_ear* column."""
    TT_stop_axis_series = pd.Series(TT_stop_axis)
    id_ear = pd.Series(1.0, index=TT_stop_axis_series.index)
    id_ear[TT_stop_axis_series.dropna().index] = np.nan
    return id_ear.tolist()


def _gen_TT_del_axis_list(TT_stop_axis_series, delais_TT_stop_del_axis):
    """Generate the *TT_del_axis* column."""
    return TT_stop_axis_series + delais_TT_stop_del_axis


def _gen_HS_final_series(axeT_, dynT_):
    """Generate the *HS_final* column."""
    HS_final_series = pd.Series(index=axeT_.index)
    dynT_grouped = dynT_.groupby(["id_axis", "N_phytomer_potential"])
    for axeT_key, axeT_group in axeT_.groupby(["id_axis", "N_phytomer_potential"]):
        dynT_group = dynT_grouped.get_group(axeT_key)
        current_a_cohort = dynT_group["a_cohort"][dynT_group.first_valid_index()]
        current_TT_hs_0 = dynT_group["TT_hs_0"][dynT_group.first_valid_index()]
        HS_final_series[axeT_group.index] = current_a_cohort * (
            axeT_group["TT_stop_axis"][axeT_group.index] - current_TT_hs_0
        )
    index_to_modify = HS_final_series[
        HS_final_series > axeT_["N_phytomer_potential"]
    ].index
    HS_final_series[index_to_modify] = axeT_["N_phytomer_potential"][index_to_modify]
    HS_final_series.fillna(axeT_["N_phytomer_potential"], inplace=True)
    HS_final_series = HS_final_series.clip(lower=0.0)
    HS_final_series = HS_final_series[HS_final_series != 0.0]
    return HS_final_series


def _remove_axes_without_leaf(axeT_, index_to_keep):
    """Remove the axes which do not have any leaf."""
    axeT_ = axeT_.loc[index_to_keep]
    axeT_.index = list(range(len(axeT_)))
    return axeT_


def _create_tilleringT(
    dynT_,
    phenT_first,
    number_of_axes,
    plants_number,
    plants_density,
    ears_density,
    TT_regression_start,
    TT_regression_end,
):
    """
    Create the :ref:`tilleringT <tilleringT>` dataframe.
    """
    dynT_most_frequent_MS = dynT_.loc[dynT_.first_valid_index()]
    id_cohort_most_frequent_MS = str(dynT_most_frequent_MS["id_cohort"])
    N_phytomer_potential_most_frequent_MS = str(
        dynT_most_frequent_MS["N_phytomer_potential"]
    ).zfill(
        2
    )  # we use N_phytomer_potential because N_phytomer_potential == N_phytomer for the most frequent MS
    id_phen_most_frequent_MS = int(
        "".join(
            [id_cohort_most_frequent_MS, N_phytomer_potential_most_frequent_MS, "1"]
        )
    )
    TT_start = phenT_first["TT_em_phytomer"][
        phenT_first[phenT_first["id_phen"] == id_phen_most_frequent_MS].index[0]
    ]

    axes_density = number_of_axes / float(plants_number) * plants_density
    return pd.DataFrame(
        {
            "TT": [TT_start, TT_regression_start, TT_regression_end],
            "axes_density": [plants_density, axes_density, ears_density],
        },
        columns=["TT", "axes_density"],
    )


def _create_cardinalityT(
    theoretical_cohort_cardinalities,
    theoretical_axis_cardinalities,
    simulated_cohorts_axes,
):
    """
    Create the :ref:`cardinalityT <cardinalityT>` dataframe.
    """
    simulated_cohort_cardinalities = (
        simulated_cohorts_axes["id_cohort"].value_counts().to_dict()
    )
    simulated_axis_cardinalities = (
        simulated_cohorts_axes.groupby(["id_cohort", "id_axis"]).size().to_dict()
    )
    cardinalityT = pd.DataFrame(
        index=list(range(len(theoretical_axis_cardinalities))),
        columns=[
            "id_cohort",
            "id_axis",
            "theoretical_cohort_cardinality",
            "simulated_cohort_cardinality",
            "theoretical_axis_cardinality",
            "simulated_axis_cardinality",
        ],
    )
    idx = 0
    for (
        id_cohort,
        id_axis,
    ), theoretical_axis_cardinality in theoretical_axis_cardinalities.items():
        cardinalityT.loc[idx,"id_cohort"] = id_cohort
        cardinalityT.loc[idx,"id_axis"] = id_axis
        cardinalityT.loc[idx,"theoretical_cohort_cardinality"] = (
            theoretical_cohort_cardinalities[id_cohort]
        )
        cardinalityT.loc[idx,"theoretical_axis_cardinality"] = theoretical_axis_cardinality
        if id_cohort in simulated_cohort_cardinalities:
            cardinalityT.loc[idx,"simulated_cohort_cardinality"] = (
                simulated_cohort_cardinalities[id_cohort]
            )
        else:
            cardinalityT.loc[idx,"simulated_cohort_cardinality"] = 0
        if (id_cohort, id_axis) in simulated_axis_cardinalities:
            cardinalityT.loc[idx,"simulated_axis_cardinality"] = (
                simulated_axis_cardinalities[(id_cohort, id_axis)]
            )
        else:
            cardinalityT.loc[idx,"simulated_axis_cardinality"] = 0
        idx += 1
    cardinalityT[
        [
            "theoretical_cohort_cardinality",
            "theoretical_axis_cardinality",
            "simulated_cohort_cardinality",
            "simulated_axis_cardinality",
        ]
    ] = cardinalityT[
        [
            "theoretical_cohort_cardinality",
            "theoretical_axis_cardinality",
            "simulated_cohort_cardinality",
            "simulated_axis_cardinality",
        ]
    ].astype(float)
    cardinalityT.loc[:,
        ["id_cohort", "simulated_cohort_cardinality", "simulated_axis_cardinality"]
    ] = cardinalityT[
        ["id_cohort", "simulated_cohort_cardinality", "simulated_axis_cardinality"]
    ].astype(int)
    cardinalityT.sort_values(["id_cohort", "id_axis"], inplace=True)
    cardinalityT.index = list(range(len(cardinalityT)))
    return cardinalityT


class _CreateDimTTmp:
    """
    Create the *dimT_tmp* dataframe.
    Compute the following columns: *id_axis*, *N_phytomer_potential*, *index_phytomer*.
    """

    def __init__(self):
        self.dimT_tmp = None

    def __call__(self, axeT_tmp, force=True):
        if force or self.dimT_tmp is None:
            id_axis_list = []
            N_phytomer_potential_list = []
            index_phytomer_list = []
            for (id_axis, N_phytomer_potential), axeT_tmp_group in axeT_tmp.groupby(
                ["id_axis", "N_phytomer_potential"]
            ):
                id_axis_list.extend(np.repeat(id_axis, N_phytomer_potential))
                N_phytomer_potential_list.extend(
                    np.repeat(N_phytomer_potential, N_phytomer_potential)
                )
                index_phytomer_list.extend(
                    list(range(1, int(N_phytomer_potential) + 1))
                )

            self.dimT_tmp = pd.DataFrame(
                index=list(range(len(id_axis_list))),
                columns=[
                    "id_axis",
                    "N_phytomer_potential",
                    "index_phytomer",
                    "L_blade",
                    "W_blade",
                    "L_sheath",
                    "W_sheath",
                    "L_internode",
                    "W_internode",
                ],
                dtype=float,
            )
            self.dimT_tmp["id_axis"] = id_axis_list
            self.dimT_tmp["N_phytomer_potential"] = N_phytomer_potential_list
            self.dimT_tmp["index_phytomer"] = index_phytomer_list
        return self.dimT_tmp


_create_dimT_tmp = _CreateDimTTmp()


class _CreateDimT:
    """
    Create the :ref:`dimT <dimT>` dataframe filling the *dimT_tmp* dataframe.
    """

    def __init__(self):
        self.dimT_ = None

    def __call__(
        self, axeT_, dimT_tmp, dynT_, decimal_elongated_internode_number, force=True
    ):
        if force or self.dimT_ is None:
            if dimT_tmp["id_axis"].count() != dimT_tmp["id_axis"].size:
                raise tools.InputError("dimT_tmp['id_axis'] contains NA values")
            if (
                dimT_tmp["N_phytomer_potential"].count()
                != dimT_tmp["N_phytomer_potential"].size
            ):
                raise tools.InputError(
                    "dimT_tmp['N_phytomer_potential'] contains NA values"
                )
            if dimT_tmp["index_phytomer"].count() != dimT_tmp["index_phytomer"].size:
                raise tools.InputError("dimT_tmp['index_phytomer'] contains NA values")

            dimT_tmp_grouped = dimT_tmp.groupby(["id_axis", "N_phytomer_potential"])
            dimT_tmp_group = dimT_tmp_grouped.get_group(
                (dynT_["id_axis"][0], dynT_["N_phytomer_potential"][0])
            )
            dimT_tmp_group_without_na = dimT_tmp_group.dropna()
            if len(dimT_tmp_group_without_na) != len(dimT_tmp_group):
                raise tools.InputError(
                    "dimT_tmp does not contain the dimensions of the most frequent MS"
                )

            self.dimT_ = _init_dimT(axeT_, dimT_tmp, dynT_)

            MS_dynT = dynT_[dynT_["id_axis"] == "MS"]
            idxmax = MS_dynT["cardinality"].idxmax()
            MS_id_cohort = MS_dynT["id_cohort"][idxmax]
            MS_N_phytomer_potential = MS_dynT["N_phytomer_potential"][idxmax]
            axeT_grouped = axeT_.groupby(["id_axis", "id_cohort", "N_phytomer"])
            axeT_group = axeT_grouped.get_group(
                ("MS", MS_id_cohort, MS_N_phytomer_potential)
            )
            MS_id_dim = axeT_group["id_dim"][axeT_group.first_valid_index()]

            L_blade_is_null = self.dimT_["L_blade"].isnull()
            row_indexes_to_fit = L_blade_is_null[L_blade_is_null].index

            _gen_lengths(
                MS_id_dim,
                row_indexes_to_fit,
                self.dimT_,
                decimal_elongated_internode_number,
            )

            _gen_widths(
                MS_id_dim,
                row_indexes_to_fit,
                self.dimT_,
                decimal_elongated_internode_number,
            )

            self.dimT_.sort_values(
                ["is_ear", "id_dim"], ascending=[False, True], inplace=True
            )

            # reinitialize the index
            self.dimT_.index = list(range(self.dimT_.index.size))
            del self.dimT_["id_cohort"]
            del self.dimT_["is_ear"]
        return self.dimT_


_create_dimT = _CreateDimT()


def _init_dimT(axeT_, dimT_tmp, dynT_):
    """Initialize dimT."""
    dimT_ = pd.DataFrame(
        columns=[
            "id_dim",
            "id_cohort",
            "index_phytomer",
            "index_relative_to_MS_phytomer",
            "L_blade",
            "W_blade",
            "L_sheath",
            "W_sheath",
            "L_internode",
            "W_internode",
            "is_ear",
        ]
    )
    dimT_tmp_grouped = dimT_tmp.groupby(["id_axis", "N_phytomer_potential"])
    for id_dim, axeT_group in axeT_.groupby("id_dim"):
        axeT_keys = list(
            axeT_group.groupby(
                ["id_axis", "id_cohort", "N_phytomer_potential"]
            ).groups.keys()
        )
        dynT_group = dynT_.loc[
            dynT_.index.map(
                lambda idx: (
                    dynT_["id_axis"][idx],
                    dynT_["id_cohort"][idx],
                    dynT_["N_phytomer_potential"][idx],
                )
                in axeT_keys
            )
        ]
        idxmax = dynT_group[
            dynT_group["id_axis"] == dynT_group["id_axis"].max()
        ].first_valid_index()
        N_phytomer_potential = dynT_group["N_phytomer_potential"][idxmax]

        dimT_group_idx = np.arange(
            axeT_group["N_phytomer"][axeT_group.first_valid_index()]
        )

        dimT_group = pd.DataFrame(
            index=dimT_group_idx, columns=dimT_.columns, dtype=float
        )
        dimT_group["id_dim"] = id_dim
        dimT_group["index_phytomer"] = dimT_group_idx + 1

        id_cohort = axeT_keys[0][1]
        dimT_group["id_cohort"] = id_cohort
        if id_cohort == 1:  # MS
            dimT_group["index_relative_to_MS_phytomer"] = dimT_group.index_phytomer
        else:
            dimT_group["index_relative_to_MS_phytomer"] = dimT_group.index_phytomer - (
                params.SLOPE_SHIFT_MS_TO_TILLERS * id_cohort
            )

        is_ear = int(str(int(id_dim))[-1])

        if is_ear == 1:
            id_axis = dynT_group["id_axis"][idxmax]
            dimT_tmp_group = dimT_tmp_grouped.get_group((id_axis, N_phytomer_potential))
            organ_dim_list = [
                "L_blade",
                "W_blade",
                "L_sheath",
                "W_sheath",
                "L_internode",
                "W_internode",
            ]
            for organ_dim in organ_dim_list:
                dim_idx_to_get = dimT_tmp_group.index[dimT_group.index.astype(int)]
                dimT_group[organ_dim] = dimT_tmp_group[organ_dim][
                    dim_idx_to_get
                ].values.astype(float)

        dimT_group["is_ear"] = is_ear

        dimT_ = pd.concat([dimT_ if not dimT_.empty else None, dimT_group], ignore_index=True)

    # force the type of id_dim, index_phytomer and is_ear
    dimT_[["id_dim", "id_cohort", "index_phytomer", "is_ear"]] = dimT_[
        ["id_dim", "id_cohort", "index_phytomer", "is_ear"]
    ].astype(int)

    return dimT_


def _gen_lengths(
    MS_id_dim, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
):
    """Fit the lengths in-place."""
    index_phytomer_series = dimT_["index_phytomer"]
    MS_rows_indexes = dimT_[dimT_["id_dim"] == MS_id_dim].index
    MS_last_but_one_index_phytomer = index_phytomer_series[MS_rows_indexes[-1] - 1]
    MS_last_index_phytomer = index_phytomer_series[MS_rows_indexes[-1]]

    _fit_L_blade(
        index_phytomer_series,
        MS_rows_indexes,
        MS_last_but_one_index_phytomer,
        MS_last_index_phytomer,
        row_indexes_to_fit,
        dimT_,
        decimal_elongated_internode_number,
    )
    _fit_L_sheath(
        index_phytomer_series,
        MS_rows_indexes,
        MS_last_but_one_index_phytomer,
        MS_last_index_phytomer,
        row_indexes_to_fit,
        dimT_,
    )
    _fit_L_internode(
        index_phytomer_series,
        MS_rows_indexes,
        MS_last_but_one_index_phytomer,
        MS_last_index_phytomer,
        row_indexes_to_fit,
        dimT_,
    )


def _fit_L_blade(
    index_phytomer_series,
    MS_rows_indexes,
    MS_last_but_one_index_phytomer,
    MS_last_index_phytomer,
    row_indexes_to_fit,
    dimT_,
    decimal_elongated_internode_number,
):
    length = "L_blade"
    lengths_series = dimT_[length]
    MS_lengths_series = lengths_series[MS_rows_indexes]
    MS_index_phytomer_series = dimT_["index_phytomer"][MS_rows_indexes]
    MS_index_phytomer_normalized_series = (
        MS_index_phytomer_series / MS_index_phytomer_series.max()
    )
    most_frequent_MS_polynomial_coefficients_array_normalized = np.polyfit(
        MS_index_phytomer_normalized_series.values, MS_lengths_series.values, 6
    )

    MS_after_start_MS_elongation_index_phytomer_series = MS_index_phytomer_series[
        MS_index_phytomer_series >= decimal_elongated_internode_number - 1
    ]
    MS_after_start_MS_elongation_polynomial_coefficients_array = np.polyfit(
        MS_after_start_MS_elongation_index_phytomer_series.values,
        MS_lengths_series[
            MS_after_start_MS_elongation_index_phytomer_series.index
        ].values,
        4,
    )

    MS_length_at_decimal_elongated_internode_number = np.polyval(
        MS_after_start_MS_elongation_polynomial_coefficients_array,
        decimal_elongated_internode_number,
    )

    MS_last_but_one_length = MS_lengths_series[MS_lengths_series.last_valid_index() - 1]
    MS_last_length = MS_lengths_series[MS_lengths_series.last_valid_index()]
    MS_last_two_lengths_coefficients_array = np.polyfit(
        [MS_last_but_one_index_phytomer, MS_last_index_phytomer],
        [MS_last_but_one_length, MS_last_length],
        1,
    )

    lengths_multiplicative_factor = 1 + params.LENGTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]

        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            index_phytomer_normalized_series = (
                index_phytomer_series / index_phytomer_series.max()
            )
            dimT_.loc[dimT_group.index, length] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array_normalized,
                index_phytomer_normalized_series.values,
            )
        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            # after start of MS elongation: polynomial phase
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_last_index_phytomer
            ].index
            after_start_MS_elongation_index_relative_to_MS_phytomer_series = (
                index_relative_to_MS_phytomer_series[
                    index_relative_to_MS_phytomer_series
                    >= decimal_elongated_internode_number
                ]
            )
            after_start_MS_elongation_index_relative_to_MS_phytomer_series_indexes = (
                after_start_MS_elongation_index_relative_to_MS_phytomer_series.index
            )
            indexes_to_compute = indexes_to_compute.intersection(
                after_start_MS_elongation_index_relative_to_MS_phytomer_series_indexes
            )
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                MS_after_start_MS_elongation_polynomial_coefficients_array,
                after_start_MS_elongation_index_relative_to_MS_phytomer_series[
                    indexes_to_compute
                ].values,
            )

            # when index_relative_to_MS_phytomer_series > MS_last_index_phytomer: same slope as MS last two lengths
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series > MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                MS_last_two_lengths_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )

            # before start of MS elongation: linear phase
            x1 = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series.first_valid_index()
            ]
            y1 = params.TILLERS_L_BLADE_1ST
            x2 = decimal_elongated_internode_number
            y2 = MS_length_at_decimal_elongated_internode_number
            before_start_MS_elongation_polynomial_coefficient_array = np.polyfit(
                np.array([x1, x2]), np.array([y1, y2]), 1
            )
            before_start_MS_elongation_index_relative_to_MS_phytomer_series = (
                index_relative_to_MS_phytomer_series[
                    index_relative_to_MS_phytomer_series
                    <= decimal_elongated_internode_number
                ]
            )
            dimT_.loc[
                before_start_MS_elongation_index_relative_to_MS_phytomer_series.index,
                length,
            ] = np.polyval(
                before_start_MS_elongation_polynomial_coefficient_array,
                before_start_MS_elongation_index_relative_to_MS_phytomer_series.values,
            )

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, length] *= lengths_multiplicative_factor


def _fit_L_sheath(
    index_phytomer_series,
    MS_rows_indexes,
    MS_last_but_one_index_phytomer,
    MS_last_index_phytomer,
    row_indexes_to_fit,
    dimT_,
):
    length = "L_sheath"
    lengths_series = dimT_[length]
    MS_lengths_series = lengths_series[MS_rows_indexes]
    MS_index_phytomer_series = dimT_["index_phytomer"][MS_rows_indexes]

    most_frequent_MS_polynomial_coefficients_array = np.polyfit(
        MS_index_phytomer_series.values, MS_lengths_series.values, 4
    )

    MS_index_phytomer_normalized_series = (
        MS_index_phytomer_series / MS_index_phytomer_series.max()
    )
    most_frequent_MS_polynomial_coefficients_array_normalized = np.polyfit(
        MS_index_phytomer_normalized_series.values, MS_lengths_series.values, 4
    )

    MS_last_but_one_length = MS_lengths_series[MS_lengths_series.last_valid_index() - 1]
    MS_last_length = MS_lengths_series[MS_lengths_series.last_valid_index()]
    MS_last_two_lengths_coefficients_array = np.polyfit(
        [MS_last_but_one_index_phytomer, MS_last_index_phytomer],
        [MS_last_but_one_length, MS_last_length],
        1,
    )

    lengths_multiplicative_factor = 1 + params.LENGTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]

        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            index_phytomer_normalized_series = (
                index_phytomer_series / index_phytomer_series.max()
            )
            dimT_.loc[dimT_group.index, length] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array_normalized,
                index_phytomer_normalized_series.values,
            )
        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )

            # when index_relative_to_MS_phytomer_series > MS_last_index_phytomer: same slope as MS last two lengths
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series > MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                MS_last_two_lengths_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, length] *= lengths_multiplicative_factor


def _fit_L_internode(
    index_phytomer_series,
    MS_rows_indexes,
    MS_last_but_one_index_phytomer,
    MS_last_index_phytomer,
    row_indexes_to_fit,
    dimT_,
):
    length = "L_internode"
    lengths_series = dimT_[length]
    MS_lengths_series = lengths_series[MS_rows_indexes]
    MS_non_null_lengths_rows_indexes = MS_lengths_series[MS_lengths_series != 0.0].index
    MS_index_phytomer_series = dimT_.loc[
        MS_non_null_lengths_rows_indexes, "index_phytomer"
    ]
    MS_first_non_null_index_phytomer = MS_index_phytomer_series[
        MS_index_phytomer_series.first_valid_index()
    ]

    most_frequent_MS_polynomial_coefficients_array = np.polyfit(
        MS_index_phytomer_series.values,
        MS_lengths_series[MS_non_null_lengths_rows_indexes].values,
        4,
    )

    MS_index_phytomer_normalized_series = (
        MS_index_phytomer_series / MS_index_phytomer_series.max()
    )
    most_frequent_MS_polynomial_coefficients_array_normalized = np.polyfit(
        MS_index_phytomer_normalized_series.values,
        MS_lengths_series[MS_non_null_lengths_rows_indexes].values,
        4,
    )

    MS_last_but_one_length = MS_lengths_series[MS_non_null_lengths_rows_indexes[-1] - 1]
    MS_last_length = MS_lengths_series[MS_non_null_lengths_rows_indexes[-1]]
    MS_last_two_lengths_coefficients_array = np.polyfit(
        [MS_last_but_one_index_phytomer, MS_last_index_phytomer],
        [MS_last_but_one_length, MS_last_length],
        1,
    )

    lengths_multiplicative_factor = 1 + params.LENGTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]

        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            # threshold
            indexes_to_threshold = index_phytomer_series[
                index_phytomer_series <= MS_first_non_null_index_phytomer
            ].index
            dimT_.loc[indexes_to_threshold, length] = 0.0
            # compute
            indexes_to_compute = index_phytomer_series.index.difference(
                indexes_to_threshold
            )
            index_phytomer_normalized_series = (
                index_phytomer_series / index_phytomer_series.max()
            )
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array_normalized,
                index_phytomer_normalized_series[indexes_to_compute].values,
            )
        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            # threshold
            indexes_to_threshold = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_first_non_null_index_phytomer
            ].index
            dimT_.loc[indexes_to_threshold, length] = 0.0
            # compute
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_last_index_phytomer
            ].index
            indexes_to_compute = indexes_to_compute.difference(indexes_to_threshold)
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )

            # when index_relative_to_MS_phytomer_series > MS_last_index_phytomer: same slope as MS last two lengths
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series > MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_compute, length] = np.polyval(
                MS_last_two_lengths_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, length] *= lengths_multiplicative_factor


def _gen_widths(
    MS_id_dim, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
):
    """Fit the widths in-place."""
    MS_rows_indexes = dimT_[dimT_["id_dim"] == MS_id_dim].index

    _fit_W_blade(MS_rows_indexes, row_indexes_to_fit, dimT_)
    _fit_W_sheath(
        MS_rows_indexes, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
    )
    _fit_W_internode(
        MS_rows_indexes, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
    )


def _fit_W_blade(MS_rows_indexes, row_indexes_to_fit, dimT_):
    MS_index_phytomer_series = dimT_["index_phytomer"][MS_rows_indexes]

    MS_first_index_phytomer = MS_index_phytomer_series[
        MS_index_phytomer_series.first_valid_index()
    ]
    MS_last_index_phytomer = MS_index_phytomer_series[
        MS_index_phytomer_series.last_valid_index()
    ]

    width = "W_blade"
    current_width_series = dimT_[width]
    MS_width_series = current_width_series[MS_rows_indexes]
    MS_first_width = MS_width_series[MS_width_series.first_valid_index()]
    tiller_first_width = MS_first_width * params.K1
    MS_last_width = MS_width_series[MS_width_series.last_valid_index()]
    tiller_last_width = MS_last_width * params.K2

    MS_index_phytomer_normalized_series = (
        MS_index_phytomer_series / MS_index_phytomer_series.max()
    )
    most_frequent_MS_polynomial_coefficients_array_normalized = np.polyfit(
        MS_index_phytomer_normalized_series.values, MS_width_series.values, 4
    )

    widths_multiplicative_factor = 1 + params.WIDTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]

        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            index_phytomer_normalized_series = (
                index_phytomer_series / index_phytomer_series.max()
            )
            dimT_.loc[dimT_group.index, width] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array_normalized,
                index_phytomer_normalized_series.values,
            )
        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            # compute
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_last_index_phytomer
            ].index
            tiller_last_index_phytomer_to_compute = dimT_group.index_phytomer.loc[
                indexes_to_compute[-1]
            ]
            most_frequent_MS_polynomial_coefficients_array = np.polyfit(
                np.array(
                    [MS_first_index_phytomer, tiller_last_index_phytomer_to_compute]
                ),
                np.array([tiller_first_width, tiller_last_width]),
                1,
            )
            dimT_.loc[indexes_to_compute, width] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array,
                index_relative_to_MS_phytomer_series[indexes_to_compute].values,
            )
            width_offset = dimT_.loc[indexes_to_compute[0], width] - MS_first_width
            dimT_.loc[indexes_to_compute, width] -= width_offset

            # ceiling
            indexes_to_ceil = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series > MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_ceil, width] = tiller_last_width

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, width] *= widths_multiplicative_factor


def _fit_W_sheath(
    MS_rows_indexes, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
):
    MS_index_phytomer_series = dimT_["index_phytomer"][MS_rows_indexes]

    MS_first_index_phytomer = MS_index_phytomer_series[
        MS_index_phytomer_series.first_valid_index()
    ]
    MS_last_index_phytomer = MS_index_phytomer_series[
        MS_index_phytomer_series.last_valid_index()
    ]

    width = "W_sheath"
    current_width_series = dimT_[width]
    MS_width_series = current_width_series[MS_rows_indexes]
    MS_first_width = MS_width_series[MS_width_series.first_valid_index()]

    # after start of MS elongation
    MS_after_start_MS_elongation_index_phytomer_series = MS_index_phytomer_series[
        MS_index_phytomer_series >= decimal_elongated_internode_number
    ]
    MS_after_start_MS_elongation_width_series = MS_width_series[
        MS_after_start_MS_elongation_index_phytomer_series.index
    ]
    mean_of_MS_width_after_start_MS_elongation = (
        MS_after_start_MS_elongation_width_series.mean()
    )

    # before start of MS elongation
    MS_before_start_MS_elongation_index_phytomer_series = MS_index_phytomer_series[
        MS_index_phytomer_series < decimal_elongated_internode_number
    ]
    MS_last_index_phytomer_before_start_MS_elongation = MS_index_phytomer_series[
        MS_before_start_MS_elongation_index_phytomer_series.last_valid_index()
    ]
    MS_before_start_MS_elongation_width_series = MS_width_series[
        MS_before_start_MS_elongation_index_phytomer_series.index
    ]
    most_frequent_MS_polynomial_coefficients_array_before_start_MS_elongation = (
        np.polyfit(
            np.array(
                [
                    MS_first_index_phytomer,
                    MS_last_index_phytomer_before_start_MS_elongation,
                ]
            ),
            np.array([MS_first_width, mean_of_MS_width_after_start_MS_elongation]),
            1,
        )
    )

    MS_width_at_decimal_elongated_internode_number = np.polyval(
        most_frequent_MS_polynomial_coefficients_array_before_start_MS_elongation,
        decimal_elongated_internode_number,
    )

    widths_multiplicative_factor = 1 + params.WIDTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]

        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            # after start of MS elongation
            after_start_MS_elongation_index_phytomer_series = index_phytomer_series[
                index_phytomer_series >= decimal_elongated_internode_number
            ]
            dimT_.loc[after_start_MS_elongation_index_phytomer_series.index, width] = (
                mean_of_MS_width_after_start_MS_elongation
            )
            # before
            before_start_MS_elongation_index_phytomer_series = index_phytomer_series[
                index_phytomer_series < decimal_elongated_internode_number
            ]
            dimT_.loc[before_start_MS_elongation_index_phytomer_series.index, width] = (
                MS_before_start_MS_elongation_width_series.values
            )
        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            # compute
            # before start of MS elongation
            before_start_MS_elongation_index_relative_to_MS_phytomer_series = (
                index_relative_to_MS_phytomer_series[
                    index_relative_to_MS_phytomer_series
                    < decimal_elongated_internode_number
                ]
            )
            dimT_.loc[
                before_start_MS_elongation_index_relative_to_MS_phytomer_series.index,
                width,
            ] = np.polyval(
                most_frequent_MS_polynomial_coefficients_array_before_start_MS_elongation,
                before_start_MS_elongation_index_relative_to_MS_phytomer_series.values,
            )

            # after start of MS elongation
            indexes_to_compute = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series <= MS_last_index_phytomer
            ].index
            after_start_MS_elongation_index_relative_to_MS_phytomer_series = (
                index_relative_to_MS_phytomer_series[
                    index_relative_to_MS_phytomer_series
                    >= decimal_elongated_internode_number
                ]
            )
            after_start_MS_elongation_index_relative_to_MS_phytomer_series_indexes = (
                after_start_MS_elongation_index_relative_to_MS_phytomer_series.index
            )
            indexes_to_compute = indexes_to_compute.intersection(
                after_start_MS_elongation_index_relative_to_MS_phytomer_series_indexes
            )
            dimT_.loc[indexes_to_compute, width] = (
                MS_width_at_decimal_elongated_internode_number
            )

            # ceiling
            indexes_to_ceil = index_relative_to_MS_phytomer_series[
                index_relative_to_MS_phytomer_series > MS_last_index_phytomer
            ].index
            dimT_.loc[indexes_to_ceil, width] = (
                mean_of_MS_width_after_start_MS_elongation
            )

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, width] *= widths_multiplicative_factor


def _fit_W_internode(
    MS_rows_indexes, row_indexes_to_fit, dimT_, decimal_elongated_internode_number
):
    width = "W_internode"

    MS_width_series = dimT_[width][MS_rows_indexes]
    MS_last_3_widths_mean = MS_width_series[-3:].mean()

    MS_N_phytomer = len(MS_rows_indexes)
    decimal_elongated_internode_number_normalized = (
        decimal_elongated_internode_number / MS_N_phytomer
    )

    most_frequent_MS_polynomial_coefficients_array = np.array(
        [
            params.W_INTERNODE_POLYNOMIAL["coefficients"]["a2"],
            params.W_INTERNODE_POLYNOMIAL["coefficients"]["a1"],
            params.W_INTERNODE_POLYNOMIAL["coefficients"]["a0"],
        ]
    )

    polynomial_first_index_relative_to_MS_phytomer_normalized = (
        params.W_INTERNODE_POLYNOMIAL["first_point"][
            "index_relative_to_MS_phytomer_normalized"
        ]
    )
    polynomial_first_W_internode_normalized = params.W_INTERNODE_POLYNOMIAL[
        "first_point"
    ]["W_internode_normalized"]

    widths_multiplicative_factor = 1 + params.WIDTHS_REDUCTION_FACTOR

    for id_dim, dimT_group in dimT_.loc[row_indexes_to_fit].groupby(by="id_dim"):
        id_cohort = dimT_group.id_cohort[dimT_group.first_valid_index()]
        if id_cohort == 1:  # MS
            index_phytomer_series = dimT_group.index_phytomer
            index_relative_to_MS_phytomer_normalized = (
                index_phytomer_series / MS_N_phytomer
            )

            # before start of MS elongation: threshold
            indexes_to_threshold = index_relative_to_MS_phytomer_normalized[
                index_relative_to_MS_phytomer_normalized
                <= decimal_elongated_internode_number_normalized
            ].index
            dimT_.loc[indexes_to_threshold, width] = 0.0

            # after start of MS elongation
            # constant phase
            indexes_to_set = index_relative_to_MS_phytomer_normalized[
                (
                    index_relative_to_MS_phytomer_normalized
                    > decimal_elongated_internode_number_normalized
                )
                & (
                    index_relative_to_MS_phytomer_normalized
                    < polynomial_first_index_relative_to_MS_phytomer_normalized
                )
            ].index
            dimT_.loc[indexes_to_set, width] = (
                polynomial_first_W_internode_normalized * MS_last_3_widths_mean
            )
            # polynomial phase
            indexes_to_compute = index_relative_to_MS_phytomer_normalized[
                (
                    index_relative_to_MS_phytomer_normalized
                    >= polynomial_first_index_relative_to_MS_phytomer_normalized
                )
                & (index_relative_to_MS_phytomer_normalized <= 1.0)
            ].index
            dimT_.loc[indexes_to_compute, width] = (
                np.polyval(
                    most_frequent_MS_polynomial_coefficients_array,
                    index_relative_to_MS_phytomer_normalized[indexes_to_compute].values,
                )
                * MS_last_3_widths_mean
            )
            # ceiling
            indexes_to_ceil = index_relative_to_MS_phytomer_normalized[
                index_relative_to_MS_phytomer_normalized > 1.0
            ].index
            dimT_.loc[indexes_to_ceil, width] = MS_last_3_widths_mean

        else:  # tiller
            index_relative_to_MS_phytomer_series = (
                dimT_group.index_relative_to_MS_phytomer
            )
            index_relative_to_MS_phytomer_normalized = (
                index_relative_to_MS_phytomer_series / MS_N_phytomer
            )

            # before start of MS elongation: threshold
            indexes_to_threshold = index_relative_to_MS_phytomer_normalized[
                index_relative_to_MS_phytomer_normalized
                <= decimal_elongated_internode_number_normalized
            ].index
            dimT_.loc[indexes_to_threshold, width] = 0.0

            # after start of MS elongation
            # constant phase
            indexes_to_set = index_relative_to_MS_phytomer_normalized[
                (
                    index_relative_to_MS_phytomer_normalized
                    > decimal_elongated_internode_number_normalized
                )
                & (
                    index_relative_to_MS_phytomer_normalized
                    < polynomial_first_index_relative_to_MS_phytomer_normalized
                )
            ].index
            dimT_.loc[indexes_to_set, width] = (
                polynomial_first_W_internode_normalized * MS_last_3_widths_mean
            )
            # polynomial phase
            indexes_to_compute = index_relative_to_MS_phytomer_normalized[
                (
                    index_relative_to_MS_phytomer_normalized
                    >= polynomial_first_index_relative_to_MS_phytomer_normalized
                )
                & (index_relative_to_MS_phytomer_normalized <= 1.0)
            ].index
            dimT_.loc[indexes_to_compute, width] = (
                np.polyval(
                    most_frequent_MS_polynomial_coefficients_array,
                    index_relative_to_MS_phytomer_normalized[indexes_to_compute].values,
                )
                * MS_last_3_widths_mean
            )
            # ceiling
            indexes_to_ceil = index_relative_to_MS_phytomer_normalized[
                index_relative_to_MS_phytomer_normalized > 1.0
            ].index
            dimT_.loc[indexes_to_ceil, width] = MS_last_3_widths_mean

            # reduction of regressive tillers
            is_ear = dimT_group.is_ear[dimT_group.first_valid_index()]
            if is_ear == 0:  # regressive
                # apply reduction factor
                dimT_.loc[dimT_group.index, width] *= widths_multiplicative_factor


def _create_dynT_tmp(axeT_tmp):
    """
    Create the *dynT_tmp* dataframe.
    """
    groups = axeT_tmp.groupby(["id_axis", "id_cohort", "N_phytomer_potential"]).groups
    keys_array = np.array(list(groups.keys()))
    cardinalities = pd.DataFrame(
        np.array(list(groups.values()), dtype=object)
    ).map(np.size)
    # initialize the values of the other columns to NaN
    dynT_tmp = pd.DataFrame(
        index=list(range(len(groups))),
        columns=[
            "id_axis",
            "id_cohort",
            "cardinality",
            "N_phytomer_potential",
            "a_cohort",
            "TT_hs_0",
            "TT_hs_break",
            "TT_flag_ligulation",
            "dTT_MS_cohort",
            "n0",
            "n1",
            "n2",
            "t0",
            "t1",
            "hs_t1",
            "a",
            "c",
            "RMSE_gl",
        ],
        dtype=float,
    )

    # set the columns 'id_axis', 'id_cohort', 'cardinality' and 'N_phytomer_potential'
    dynT_tmp["id_axis"] = keys_array[:, 0]
    dynT_tmp["id_cohort"] = keys_array[:, 1].astype(float).astype(int)
    dynT_tmp["cardinality"] = cardinalities
    dynT_tmp["N_phytomer_potential"] = keys_array[:, 2].astype(float).astype(int)

    # nested sort of dynT_tmp: first by 'id_axis' in ascending order, second
    # by 'cardinality' in descending order.
    dynT_tmp.sort_values(["id_axis", "cardinality"], ascending=[1, 0], inplace=True)
    # reinitialize the index
    dynT_tmp.index = list(range(dynT_tmp.index.size))
    return dynT_tmp


def _create_dynT(
    dynT_tmp,
    GL_number,
    leaf_number_delay_MS_cohort=params.MS_HS_AT_TILLER_EMERGENCE,
    TT_t1_user=None,
):
    """
    Create the :ref:`dynT <dynT>` dataframe.
    """
    # in 'dynT_tmp', check that the columns 'id_axis', 'cardinality' and
    # 'N_phytomer_potential' are non-NA.
    if dynT_tmp["id_axis"].count() != dynT_tmp["id_axis"].size:
        raise tools.InputError("dynT_tmp['id_axis'] contains NA values")
    if dynT_tmp["cardinality"].count() != dynT_tmp["cardinality"].size:
        raise tools.InputError("dynT_tmp['cardinality'] contains NA values")
    if (
        dynT_tmp["N_phytomer_potential"].count()
        != dynT_tmp["N_phytomer_potential"].size
    ):
        raise tools.InputError("dynT_tmp['N_phytomer_potential'] contains NA values")
    # get all main stem rows
    MS = dynT_tmp[dynT_tmp["id_axis"] == "MS"]

    # get the row of the most frequent main stem
    most_frequent_MS = MS.loc[0:0]
    # for this row, fill the columns referring to the dynamic of the green leaves
    most_frequent_MS = _gen_most_frequent_MS_GL_dynamic(
        most_frequent_MS, GL_number, TT_t1_user
    )

    # get the rows of all main stems except the most frequent one
    other_MS = MS.loc[1:]
    # for these rows, fill the columns referring to the dynamic of Haun Stage
    other_MS = _gen_other_MS_HS_dynamic(most_frequent_MS, other_MS)
    # and fill the columns referring to the dynamic of the green leaves
    other_MS = _gen_other_MS_GL_dynamic(most_frequent_MS, other_MS)

    # get the rows corresponding to the tiller axes
    tiller_axes = dynT_tmp[dynT_tmp["id_axis"] != "MS"]
    # degenerated case of one plant without tillers
    if len(tiller_axes) <= 0:
        dynT_ = pd.concat([most_frequent_MS, other_MS])
    else:
        # extract the rows corresponding to the most frequent tiller axes
        grouped = tiller_axes.groupby("id_axis")
        most_frequent_tiller_axes = []
        for id_axis, group_indexes in grouped.groups.items():
            most_frequent_tiller_axes.append(tiller_axes.loc[group_indexes[0:1]])
        # concatenate these rows in one dataframe ; 'most_frequent_tiller_axes' is
        # now a pd.DataFrame (and is not a list anymore)
        most_frequent_tiller_axes = pd.concat(most_frequent_tiller_axes)
        # in this new dataframe, fill the columns referring to the dynamic of Haun Stage
        most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_HS_dynamic(
            most_frequent_MS, most_frequent_tiller_axes, leaf_number_delay_MS_cohort
        )
        # and fill the columns referring to the dynamic of the green leaves
        most_frequent_tiller_axes = _gen_most_frequent_tiller_axes_GL_dynamic(
            most_frequent_MS, most_frequent_tiller_axes
        )
        # extract the rows corresponding to all tiller axes except the most frequent ones
        other_tiller_axes = tiller_axes.drop(most_frequent_tiller_axes.index)
        other_tiller_axes = _gen_other_tiller_axes_HS_dynamic(
            most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes
        )
        # and fill the columns referring to the dynamic of the green leaves
        other_tiller_axes = _gen_other_tiller_axes_GL_dynamic(
            most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes
        )

        dynT_ = pd.concat(
            [most_frequent_MS, other_MS, most_frequent_tiller_axes, other_tiller_axes]
        )

    dynT_.sort_values(["id_axis", "cardinality"], ascending=[1, 0], inplace=True)
    # reinitialize the index
    dynT_.index = list(range(dynT_.index.size))

    return dynT_


def _gen_most_frequent_MS_GL_dynamic(most_frequent_MS, GL_number, TT_t1_user=None):
    """
    Create a copy of *most_frequent_MS*, fill this copy by calculating the
    parameters which describe the dynamic of the green leaves, and return this
    copy.
    """
    # calculation of t1
    most_frequent_MS = most_frequent_MS.copy()

    if TT_t1_user is None:
        MS_HS_for_minimum_GL = (
            most_frequent_MS["N_phytomer_potential"][0]
            - params.NUMBER_OF_ELONGATED_INTERNODES
        )
        if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
            most_frequent_MS["t1"] = (
                MS_HS_for_minimum_GL / most_frequent_MS["a_cohort"]
                + most_frequent_MS["TT_hs_0"]
            )
        else:  # bilinear mode
            HS_break_0 = most_frequent_MS["a_cohort"][0] * (
                most_frequent_MS["TT_hs_break"][0] - most_frequent_MS["TT_hs_0"][0]
            )
            a2_0 = (most_frequent_MS["N_phytomer_potential"][0] - HS_break_0) / (
                most_frequent_MS["TT_flag_ligulation"][0]
                - most_frequent_MS["TT_hs_break"][0]
            )
            if MS_HS_for_minimum_GL < HS_break_0:
                most_frequent_MS["t1"] = (
                    most_frequent_MS["TT_hs_0"]
                    + MS_HS_for_minimum_GL / most_frequent_MS["a_cohort"]
                )
            else:
                most_frequent_MS["t1"] = (
                    MS_HS_for_minimum_GL - HS_break_0
                ) / a2_0 + most_frequent_MS["TT_hs_break"]
    else:
        most_frequent_MS["t1"] = TT_t1_user
    # calculation of hs_t1
    most_frequent_MS["hs_t1"] = most_frequent_MS["a_cohort"] * (
        most_frequent_MS["t1"] - most_frequent_MS["TT_hs_0"]
    )
    # calculation of t0
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
        most_frequent_MS["t0"] = (
            most_frequent_MS["TT_hs_0"]
            + most_frequent_MS["n0"] / most_frequent_MS["a_cohort"]
        )
    else:  # bilinear mode
        HS_break = most_frequent_MS["a_cohort"] * (
            most_frequent_MS["TT_hs_break"] - most_frequent_MS["TT_hs_0"]
        )
        a2 = (most_frequent_MS["N_phytomer_potential"] - HS_break) / (
            most_frequent_MS["TT_flag_ligulation"] - most_frequent_MS["TT_hs_break"]
        )
        n0_smaller_than_HS_break_indexes = most_frequent_MS[
            most_frequent_MS["n0"] < HS_break
        ].index
        n0_greater_than_HS_break_indexes = most_frequent_MS[
            most_frequent_MS["n0"] >= HS_break
        ].index
        most_frequent_MS.loc[n0_smaller_than_HS_break_indexes, "t0"] = (
            most_frequent_MS["TT_hs_0"][n0_smaller_than_HS_break_indexes]
            + most_frequent_MS["n0"][n0_smaller_than_HS_break_indexes]
            / most_frequent_MS["a_cohort"][n0_smaller_than_HS_break_indexes]
        )
        most_frequent_MS.loc[n0_greater_than_HS_break_indexes, "t0"] = (
            most_frequent_MS["n0"][n0_greater_than_HS_break_indexes]
            - HS_break[n0_greater_than_HS_break_indexes]
        ) / a2[n0_greater_than_HS_break_indexes] + most_frequent_MS["TT_hs_break"][
            n0_greater_than_HS_break_indexes
        ]
    # calculation of c
    most_frequent_MS["c"] = -(
        (most_frequent_MS["N_phytomer_potential"] - most_frequent_MS["hs_t1"])
        - (most_frequent_MS["n2"] - most_frequent_MS["n1"])
    ) / (most_frequent_MS["TT_flag_ligulation"] - most_frequent_MS["t1"])
    # calculation of a
    TT_flag_ligulation_0 = most_frequent_MS["TT_flag_ligulation"][0]
    n2_0 = most_frequent_MS["n2"][0]
    TT = (
        np.array([TT_flag_ligulation_0] + list(GL_number.keys())) - TT_flag_ligulation_0
    )
    GL = np.array([n2_0] + list(GL_number.values()))
    fixed_coefs = [0.0, most_frequent_MS["c"][0], n2_0]
    a_starting_estimate = -4.0e-9
    a_tmp, RMSE_gl = tools.fit_poly(TT, GL, fixed_coefs, a_starting_estimate)
    most_frequent_MS.loc[0, "a"], most_frequent_MS.loc[0, "RMSE_gl"] = (
        -abs(a_tmp),
        RMSE_gl,
    )  # force 'a' to be negative
    return most_frequent_MS


def _gen_other_MS_HS_dynamic(most_frequent_MS, other_MS):
    """
    Create a copy of *other_MS*, fill this copy by calculating the
    parameters which describe the dynamic of the Haun Stage, and return this
    copy.
    """
    other_MS = other_MS.copy()
    if other_MS["TT_hs_0"].count() != other_MS["TT_hs_0"].size:
        # calculation of TT_hs_0
        other_MS["TT_hs_0"] = most_frequent_MS["TT_hs_0"][0]
        # calculation of TT_hs_break
        other_MS["TT_hs_break"] = most_frequent_MS["TT_hs_break"][0]
        # calculation of dTT_MS_cohort
        other_MS["dTT_MS_cohort"] = (
            most_frequent_MS["dTT_MS_cohort"][0]
            + (
                other_MS["N_phytomer_potential"]
                - most_frequent_MS["N_phytomer_potential"][0]
            )
            * params.FLAG_LIGULATION_DELAY
        )
        # calculation of TT_flag_ligulation
        other_MS["TT_flag_ligulation"] = most_frequent_MS["TT_flag_ligulation"][
            0
        ] + other_MS["N_phytomer_potential"] / (
            most_frequent_MS["N_phytomer_potential"][0]
            * 4
            * most_frequent_MS["a_cohort"][0]
        )
        # calculation of a_cohort
        if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
            other_MS["a_cohort"] = other_MS["N_phytomer_potential"] / (
                other_MS["TT_flag_ligulation"] - other_MS["TT_hs_0"]
            )
        else:  # bilinear mode
            other_MS["a_cohort"] = most_frequent_MS["a_cohort"][0]
    return other_MS


def _gen_other_MS_GL_dynamic(most_frequent_MS, other_MS):
    """
    Create a copy of *other_MS*, fill this copy by calculating the
    parameters which describe the dynamic of the green leaves, and return this
    copy.
    """
    other_MS = other_MS.copy()
    # calculation of n1
    if other_MS["n1"].count() != other_MS["n1"].size:
        other_MS["n1"] = (
            most_frequent_MS["n1"][0]
            * other_MS["N_phytomer_potential"]
            / most_frequent_MS["N_phytomer_potential"][0]
        )
    # calculation of t1
    other_MS["t1"] = (
        other_MS["N_phytomer_potential"] - params.NUMBER_OF_ELONGATED_INTERNODES
    ) / other_MS["a_cohort"] + other_MS["TT_hs_0"]
    # calculation of hs_t1
    other_MS["hs_t1"] = other_MS["a_cohort"] * (other_MS["t1"] - other_MS["TT_hs_0"])
    # calculation of n0
    if other_MS["n0"].count() != other_MS["n0"].size:
        other_MS["n0"] = (
            most_frequent_MS["n0"][0]
            * other_MS["N_phytomer_potential"]
            / most_frequent_MS["N_phytomer_potential"][0]
        )
    # calculation of n2
    if other_MS["n2"].count() != other_MS["n2"].size:
        other_MS["n2"] = (
            most_frequent_MS["n2"][0]
            * other_MS["N_phytomer_potential"]
            / most_frequent_MS["N_phytomer_potential"][0]
        )
    # calculation of t0
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
        other_MS["t0"] = other_MS["TT_hs_0"] + other_MS["n0"] / other_MS["a_cohort"]
    else:  # bilinear mode
        HS_break = other_MS["a_cohort"] * (
            other_MS["TT_hs_break"] - other_MS["TT_hs_0"]
        )
        a2 = (other_MS["N_phytomer_potential"] - HS_break) / (
            other_MS["TT_flag_ligulation"] - other_MS["TT_hs_break"]
        )
        n0_smaller_than_HS_break_indexes = other_MS[other_MS["n0"] < HS_break].index
        n0_greater_than_HS_break_indexes = other_MS[other_MS["n0"] >= HS_break].index
        other_MS.loc[n0_smaller_than_HS_break_indexes, "t0"] = (
            other_MS["TT_hs_0"][n0_smaller_than_HS_break_indexes]
            + other_MS["n0"][n0_smaller_than_HS_break_indexes]
            / other_MS["a_cohort"][n0_smaller_than_HS_break_indexes]
        )
        other_MS.loc[n0_greater_than_HS_break_indexes, "t0"] = (
            other_MS["n0"][n0_greater_than_HS_break_indexes]
            - HS_break[n0_greater_than_HS_break_indexes]
        ) / a2[n0_greater_than_HS_break_indexes] + other_MS["TT_hs_break"][
            n0_greater_than_HS_break_indexes
        ]
    # calculation of c
    other_MS["c"] = (
        most_frequent_MS["c"][0] * other_MS["n2"] / most_frequent_MS["n2"][0]
    )
    # calculation of a
    other_MS["a"] = (
        most_frequent_MS["a"][0] * other_MS["n2"] / most_frequent_MS["n2"][0]
    )
    # calculation of RMSE_gl
    other_MS["RMSE_gl"] = most_frequent_MS["RMSE_gl"][0]
    return other_MS


def _gen_most_frequent_tiller_axes_HS_dynamic(
    most_frequent_MS, most_frequent_tiller_axes, leaf_number_delay_MS_cohort
):
    """
    Create a copy of *most_frequent_tiller_axes*, fill this copy by calculating the
    parameters which describe the dynamic of the Haun Stage, and return this
    copy.
    """
    most_frequent_tiller_axes = most_frequent_tiller_axes.copy()
    # calculation of TT_hs_break
    most_frequent_tiller_axes["TT_hs_break"] = most_frequent_MS["TT_hs_break"][0]
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.dropna(
        subset=["TT_hs_0"]
    ).index
    try:
        nan_most_frequent_tiller_axis_indexes = (
            most_frequent_tiller_axes.index.difference(
                without_nan_most_frequent_tiller_axis_indexes
            )
        )
    except AttributeError:  # backward compatibility with pandas < 0.16
        nan_most_frequent_tiller_axis_indexes = (
            most_frequent_tiller_axes.index
            - without_nan_most_frequent_tiller_axis_indexes
        )
    if len(nan_most_frequent_tiller_axis_indexes) == 0:
        return most_frequent_tiller_axes
    # calculation of TT_hs_0
    cohorts = most_frequent_tiller_axes["id_axis"][
        nan_most_frequent_tiller_axis_indexes
    ].values
    primary_cohorts = np.apply_along_axis(
        np.vectorize(tools.get_primary_axis), 0, cohorts, [params.FIRST_CHILD_DELAY]
    )
    leaf_number_delay_MS_cohorts = np.array(
        [leaf_number_delay_MS_cohort[cohort] for cohort in primary_cohorts]
    )
    most_frequent_tiller_axes.loc[nan_most_frequent_tiller_axis_indexes, "TT_hs_0"] = (
        most_frequent_MS["TT_hs_0"][0]
        + (leaf_number_delay_MS_cohorts / most_frequent_MS["a_cohort"][0])
    )
    # dTT_MS_cohort is set by the user. Thus there is nothing to do.
    # calculation of TT_flag_ligulation
    most_frequent_tiller_axes.loc[
        nan_most_frequent_tiller_axis_indexes, "TT_flag_ligulation"
    ] = (
        most_frequent_tiller_axes["dTT_MS_cohort"][
            nan_most_frequent_tiller_axis_indexes
        ]
        + most_frequent_MS["TT_flag_ligulation"][0]
    )
    # calculation of a_cohort
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
        most_frequent_tiller_axes.loc[
            nan_most_frequent_tiller_axis_indexes, "a_cohort"
        ] = most_frequent_tiller_axes["N_phytomer_potential"][
            nan_most_frequent_tiller_axis_indexes
        ] / (
            most_frequent_tiller_axes["TT_flag_ligulation"][
                nan_most_frequent_tiller_axis_indexes
            ]
            - most_frequent_tiller_axes["TT_hs_0"][
                nan_most_frequent_tiller_axis_indexes
            ]
        )
    else:  # bilinear mode
        most_frequent_tiller_axes.loc[
            nan_most_frequent_tiller_axis_indexes, "a_cohort"
        ] = most_frequent_MS["a_cohort"][0]
    return most_frequent_tiller_axes


def _gen_most_frequent_tiller_axes_GL_dynamic(
    most_frequent_MS, most_frequent_tiller_axes
):
    """
    Create a copy of *most_frequent_tiller_axes*, fill this copy by calculating the
    parameters which describe the dynamic of the green leaves, and return this
    copy.
    """
    most_frequent_tiller_axes = most_frequent_tiller_axes.copy()
    without_nan_most_frequent_tiller_axis_indexes = most_frequent_tiller_axes.dropna(
        subset=["n1"]
    ).index
    try:
        nan_most_frequent_tiller_axis_indexes = (
            most_frequent_tiller_axes.index.difference(
                without_nan_most_frequent_tiller_axis_indexes
            )
        )
    except AttributeError:  # backward compatibility with pandas < 0.16
        nan_most_frequent_tiller_axis_indexes = (
            most_frequent_tiller_axes.index
            - without_nan_most_frequent_tiller_axis_indexes
        )
    # calculation of n1
    most_frequent_tiller_axes.loc[nan_most_frequent_tiller_axis_indexes, "n1"] = (
        most_frequent_MS["n1"][0]
    )
    # calculation of t1
    most_frequent_tiller_axes["t1"] = (
        most_frequent_MS["t1"][0] + most_frequent_tiller_axes["dTT_MS_cohort"]
    )
    # calculation of hs_t1
    most_frequent_tiller_axes["hs_t1"] = most_frequent_tiller_axes["a_cohort"] * (
        most_frequent_tiller_axes["t1"] - most_frequent_tiller_axes["TT_hs_0"]
    )
    # calculation of n0
    most_frequent_tiller_axes.loc[nan_most_frequent_tiller_axis_indexes, "n2"] = (
        most_frequent_MS["n2"][0] * params.N2_MS_DIV_N2_COHORT
    )
    # calculation of t0 and n0
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
        GL_MS_line = (
            (most_frequent_MS["t0"][0], most_frequent_MS["n0"][0]),
            (most_frequent_MS["t1"][0], most_frequent_MS["n1"][0]),
        )
        for row_index, row in most_frequent_tiller_axes.iterrows():
            HS_tiller_line = (
                (row["TT_hs_0"], 0.0),
                (
                    row["TT_flag_ligulation"],
                    row["a_cohort"] * (row["TT_flag_ligulation"] - row["TT_hs_0"]),
                ),
            )
            t0_tiller, n0_tiller = tools.find_lines_intersection(
                GL_MS_line, HS_tiller_line
            )
            if t0_tiller >= most_frequent_MS["t1"][0]:
                t0_tiller = n0_tiller = np.nan
                t1_tiller = most_frequent_MS["n1"][0] / row["a_cohort"] + row["TT_hs_0"]
            else:
                t1_tiller = row["t1"]
            most_frequent_tiller_axes.loc[row_index, ["t0", "n0", "t1"]] = (
                t0_tiller,
                n0_tiller,
                t1_tiller,
            )
    else:  # bilinear mode
        GL_MS_line = (
            (most_frequent_MS["t0"][0], most_frequent_MS["n0"][0]),
            (most_frequent_MS["t1"][0], most_frequent_MS["n1"][0]),
        )
        HS_break = most_frequent_tiller_axes["a_cohort"] * (
            most_frequent_tiller_axes["TT_hs_break"]
            - most_frequent_tiller_axes["TT_hs_0"]
        )
        a2 = (most_frequent_tiller_axes["N_phytomer_potential"] - HS_break) / (
            most_frequent_tiller_axes["TT_flag_ligulation"]
            - most_frequent_tiller_axes["TT_hs_break"]
        )
        for row_index, row in most_frequent_tiller_axes.iterrows():
            HS_tiller_line_before_TT_hs_break = (
                (row["TT_hs_0"], 0.0),
                (
                    row["TT_hs_break"],
                    row["a_cohort"] * (row["TT_hs_break"] - row["TT_hs_0"]),
                ),
            )
            t0_tiller_before_TT_hs_break, n0_tiller_before_TT_hs_break = (
                tools.find_lines_intersection(
                    GL_MS_line, HS_tiller_line_before_TT_hs_break
                )
            )
            HS_tiller_line_after_TT_hs_break = (
                (row["TT_hs_break"], HS_break[row_index]),
                (
                    row["TT_flag_ligulation"],
                    a2[row_index] * (row["TT_flag_ligulation"] - row["TT_hs_break"])
                    + HS_break[row_index],
                ),
            )
            t0_tiller_after_TT_hs_break, n0_tiller_after_TT_hs_break = (
                tools.find_lines_intersection(
                    GL_MS_line, HS_tiller_line_after_TT_hs_break
                )
            )
            if t0_tiller_before_TT_hs_break <= t0_tiller_after_TT_hs_break:
                t0_tiller = t0_tiller_after_TT_hs_break
                n0_tiller = n0_tiller_after_TT_hs_break
            else:
                t0_tiller = t0_tiller_before_TT_hs_break
                n0_tiller = n0_tiller_before_TT_hs_break
            if t0_tiller >= most_frequent_MS["t1"][0]:
                t0_tiller = n0_tiller = np.nan
                t1_tiller = (most_frequent_MS["n1"][0] - HS_break[row_index]) / a2[
                    row_index
                ] + most_frequent_MS["TT_hs_break"][0]
            else:
                t1_tiller = row["t1"]
                if most_frequent_MS["TT_hs_break"][0] > most_frequent_MS["t1"][0]:
                    n0_tiller = min(n0_tiller, most_frequent_MS["n1"][0])
            most_frequent_tiller_axes.loc[row_index, ["t0", "n0", "t1"]] = (
                t0_tiller,
                n0_tiller,
                t1_tiller,
            )

    # calculation of c
    most_frequent_tiller_axes["c"] = (
        most_frequent_MS["c"][0]
        * most_frequent_tiller_axes["n2"]
        / most_frequent_MS["n2"][0]
    )
    # calculation of a
    most_frequent_tiller_axes["a"] = (
        most_frequent_MS["a"][0]
        * most_frequent_tiller_axes["n2"]
        / most_frequent_MS["n2"][0]
    )
    # calculation of RMSE_gl
    most_frequent_tiller_axes["RMSE_gl"] = most_frequent_MS["RMSE_gl"][0]
    return most_frequent_tiller_axes


def _gen_other_tiller_axes_HS_dynamic(
    most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes
):
    """
    Create a copy of *other_tiller_axes*, fill this copy by calculating the
    parameters which describe the dynamic of the Haun Stage, and return this
    copy.
    """
    other_tiller_axes = other_tiller_axes.copy()
    # calculation of TT_hs_break
    other_tiller_axes["TT_hs_break"] = most_frequent_MS["TT_hs_break"][0]
    without_nan_other_tiller_axis_indexes = other_tiller_axes.dropna(
        subset=["TT_hs_0"]
    ).index
    try:
        nan_other_tiller_axis_indexes = other_tiller_axes.index.difference(
            without_nan_other_tiller_axis_indexes
        )
    except AttributeError:  # backward compatibility with pandas < 0.16
        nan_other_tiller_axis_indexes = (
            other_tiller_axes.index - without_nan_other_tiller_axis_indexes
        )
    for name, group in other_tiller_axes.loc[nan_other_tiller_axis_indexes].groupby(
        "id_axis"
    ):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes[
            most_frequent_tiller_axes["id_axis"] == name
        ].first_valid_index()
        # calculation of TT_hs_0
        other_tiller_axes.loc[group.index, "TT_hs_0"] = most_frequent_tiller_axes[
            "TT_hs_0"
        ][most_frequent_tiller_axis_idx]
        # calculation of dTT_MS_cohort
        other_tiller_axes.loc[group.index, "dTT_MS_cohort"] = (
            most_frequent_tiller_axes["dTT_MS_cohort"][most_frequent_tiller_axis_idx]
            + (
                group["N_phytomer_potential"]
                - most_frequent_tiller_axes["N_phytomer_potential"][
                    most_frequent_tiller_axis_idx
                ]
            )
            * params.FLAG_LIGULATION_DELAY
        )
        if not math.isnan(most_frequent_MS["TT_hs_break"][0]):
            # calculation of a_cohort in bilinear mode
            other_tiller_axes.loc[group.index, "a_cohort"] = most_frequent_tiller_axes[
                "a_cohort"
            ][most_frequent_tiller_axis_idx]
    # calculation of TT_flag_ligulation
    other_tiller_axes.loc[nan_other_tiller_axis_indexes, "TT_flag_ligulation"] = (
        other_tiller_axes["dTT_MS_cohort"][nan_other_tiller_axis_indexes]
        + most_frequent_MS["TT_flag_ligulation"][0]
    )
    # calculation of a_cohort in linear mode
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):
        other_tiller_axes.loc[nan_other_tiller_axis_indexes, "a_cohort"] = (
            other_tiller_axes["N_phytomer_potential"][nan_other_tiller_axis_indexes]
            / (
                other_tiller_axes["TT_flag_ligulation"][nan_other_tiller_axis_indexes]
                - other_tiller_axes["TT_hs_0"][nan_other_tiller_axis_indexes]
            )
        )
    return other_tiller_axes


def _gen_other_tiller_axes_GL_dynamic(
    most_frequent_MS, most_frequent_tiller_axes, other_tiller_axes
):
    """
    Create a copy of *other_tiller_axes*, fill this copy by calculating the
    parameters which describe the dynamic of the green leaves, and return this
    copy.
    """
    other_tiller_axes = other_tiller_axes.copy()
    without_nan_other_tiller_axis_indexes = other_tiller_axes.dropna(
        subset=["n1"]
    ).index
    try:
        nan_other_tiller_axis_indexes = other_tiller_axes.index.difference(
            without_nan_other_tiller_axis_indexes
        )
    except AttributeError:  # backward compatibility with pandas < 0.16
        nan_other_tiller_axis_indexes = (
            other_tiller_axes.index - without_nan_other_tiller_axis_indexes
        )
    for name, group in other_tiller_axes.loc[nan_other_tiller_axis_indexes].groupby(
        "id_axis"
    ):
        most_frequent_tiller_axis_idx = most_frequent_tiller_axes[
            most_frequent_tiller_axes["id_axis"] == name
        ].first_valid_index()
        # calculation ofn1
        other_tiller_axes.loc[group.index, "n1"] = most_frequent_tiller_axes["n1"][
            most_frequent_tiller_axis_idx
        ]
        # calculation of n2
        other_tiller_axes.loc[group.index, "n2"] = most_frequent_tiller_axes["n2"][
            most_frequent_tiller_axis_idx
        ]
    # calculation of t1
    other_tiller_axes["t1"] = (
        most_frequent_MS["t1"][0] + other_tiller_axes["dTT_MS_cohort"]
    )
    # calculation of hs_t1
    other_tiller_axes["hs_t1"] = other_tiller_axes["a_cohort"] * (
        other_tiller_axes["t1"] - other_tiller_axes["TT_hs_0"]
    )
    # calculation of t0 and n0
    if math.isnan(most_frequent_MS["TT_hs_break"][0]):  # linear mode
        GL_MS_line = (
            (most_frequent_MS["t0"][0], most_frequent_MS["n0"][0]),
            (most_frequent_MS["t1"][0], most_frequent_MS["n1"][0]),
        )
        for row_index, row in other_tiller_axes.iterrows():
            HS_tiller_line = (
                (row["TT_hs_0"], 0.0),
                (
                    row["TT_flag_ligulation"],
                    row["a_cohort"] * (row["TT_flag_ligulation"] - row["TT_hs_0"]),
                ),
            )
            t0_tiller, n0_tiller = tools.find_lines_intersection(
                GL_MS_line, HS_tiller_line
            )
            if t0_tiller >= most_frequent_MS["t1"][0]:
                t0_tiller = n0_tiller = np.nan
                t1_tiller = most_frequent_MS["n1"][0] / row["a_cohort"] + row["TT_hs_0"]
            else:
                t1_tiller = row["t1"]
            other_tiller_axes.loc[row_index, ["t0", "n0", "t1"]] = (
                t0_tiller,
                n0_tiller,
                t1_tiller,
            )
    else:  # bilinear mode
        GL_MS_line = (
            (most_frequent_MS["t0"][0], most_frequent_MS["n0"][0]),
            (most_frequent_MS["t1"][0], most_frequent_MS["n1"][0]),
        )
        HS_break = other_tiller_axes["a_cohort"] * (
            other_tiller_axes["TT_hs_break"] - other_tiller_axes["TT_hs_0"]
        )
        a2 = (other_tiller_axes["N_phytomer_potential"] - HS_break) / (
            other_tiller_axes["TT_flag_ligulation"] - other_tiller_axes["TT_hs_break"]
        )
        for row_index, row in other_tiller_axes.iterrows():
            HS_tiller_line_before_TT_hs_break = (
                (row["TT_hs_0"], 0.0),
                (
                    row["TT_hs_break"],
                    row["a_cohort"] * (row["TT_hs_break"] - row["TT_hs_0"]),
                ),
            )
            t0_tiller_before_TT_hs_break, n0_tiller_before_TT_hs_break = (
                tools.find_lines_intersection(
                    GL_MS_line, HS_tiller_line_before_TT_hs_break
                )
            )
            HS_tiller_line_after_TT_hs_break = (
                (row["TT_hs_break"], HS_break[row_index]),
                (
                    row["TT_flag_ligulation"],
                    a2[row_index] * (row["TT_flag_ligulation"] - row["TT_hs_break"])
                    + HS_break[row_index],
                ),
            )
            t0_tiller_after_TT_hs_break, n0_tiller_after_TT_hs_break = (
                tools.find_lines_intersection(
                    GL_MS_line, HS_tiller_line_after_TT_hs_break
                )
            )
            if t0_tiller_before_TT_hs_break <= t0_tiller_after_TT_hs_break:
                t0_tiller = t0_tiller_after_TT_hs_break
                n0_tiller = n0_tiller_after_TT_hs_break
            else:
                t0_tiller = t0_tiller_before_TT_hs_break
                n0_tiller = n0_tiller_before_TT_hs_break

            if t0_tiller >= most_frequent_MS["t1"][0]:
                t0_tiller = n0_tiller = np.nan
                t1_tiller = (most_frequent_MS["n1"][0] - HS_break[row_index]) / a2[
                    row_index
                ] + most_frequent_MS["TT_hs_break"][0]
            else:
                t1_tiller = row["t1"]
                if most_frequent_MS["TT_hs_break"][0] > most_frequent_MS["t1"][0]:
                    n0_tiller = min(n0_tiller, most_frequent_MS["n1"][0])
            other_tiller_axes.loc[row_index, ["t0", "n0", "t1"]] = (
                t0_tiller,
                n0_tiller,
                t1_tiller,
            )
    # calculation of c
    other_tiller_axes["c"] = (
        most_frequent_MS["c"][0] * other_tiller_axes["n2"] / most_frequent_MS["n2"][0]
    )
    # calculation of a
    other_tiller_axes["a"] = (
        most_frequent_MS["a"][0] * other_tiller_axes["n2"] / most_frequent_MS["n2"][0]
    )
    # calculation of RMSE_gl
    other_tiller_axes["RMSE_gl"] = most_frequent_MS["RMSE_gl"][0]
    return other_tiller_axes


def _calculate_decimal_elongated_internode_number(dimT_tmp, dynT_tmp):
    """
    Calculate the number of elongated internodes.
    """
    MS_most_frequent_axis_dynT_tmp = dynT_tmp.loc[dynT_tmp.first_valid_index()]
    id_axis = MS_most_frequent_axis_dynT_tmp["id_axis"]
    N_phytomer_potential = MS_most_frequent_axis_dynT_tmp["N_phytomer_potential"]
    grouped = dimT_tmp.groupby(["id_axis", "N_phytomer_potential"])
    MS_most_frequent_axis_indexes = grouped.groups[(id_axis, N_phytomer_potential)]
    # get the lengths of the internodes which belong to each phytomer of the most
    # frequent axis of the MS
    MS_most_frequent_axis_L_internode = dimT_tmp["L_internode"][
        MS_most_frequent_axis_indexes
    ]
    # keep only the non-zero lengths
    MS_most_frequent_axis_L_internode = MS_most_frequent_axis_L_internode[
        MS_most_frequent_axis_L_internode != 0.0
    ]
    # get the indexes of the phytomers of the MS most frequent axis, for which L_internode is non-null
    MS_most_frequent_axis_index_phytomer = dimT_tmp["index_phytomer"][
        MS_most_frequent_axis_L_internode.index
    ]
    # Fit a polynomial of degree 2 to points (MS_most_frequent_axis_L_internode,
    # MS_most_frequent_axis_index_phytomer), and get the coefficient of degree 0.
    return np.polyfit(
        MS_most_frequent_axis_L_internode, MS_most_frequent_axis_index_phytomer, 2
    )[2]


class _CreatePhenTTmp:
    """
    Create the *phenT_tmp* dataframe.
    Compute all the columns, but the column 'TT_del_phytomer' is temporary and will
    be recalculated in _create_phenT_abs using dimT_.
    In dataframe *axeT_*, the following columns must be fulfilled: 'id_phen', 'id_cohort', 'N_phytomer'.
    """

    def __init__(self):
        self.phenT_tmp = None

    def __call__(self, axeT_, dynT_, decimal_elongated_internode_number, force=True):
        if force or self.phenT_tmp is None:
            columns_without_nan = [
                col for col in dynT_.columns if col not in ["TT_hs_break", "n0", "t0"]
            ]
            dynT_without_nan = dynT_.loc[:, columns_without_nan]
            if not (
                dynT_without_nan.count().max()
                == dynT_without_nan.count().min()
                == dynT_without_nan.index.size
            ):
                raise tools.InputError("dynT contains unexpected NA values")

            id_phen_list = []
            index_phytomer_list = []
            for (id_phen, N_phytomer), axeT_group in axeT_.groupby(
                ["id_phen", "N_phytomer"]
            ):
                id_phen_list.extend(np.repeat(id_phen, int(N_phytomer) + 1))
                index_phytomer_list.extend(list(range(int(N_phytomer) + 1)))

            self.phenT_tmp = pd.DataFrame(
                index=list(range(len(id_phen_list))),
                columns=[
                    "id_phen",
                    "index_phytomer",
                    "TT_em_phytomer",
                    "TT_col_phytomer",
                    "TT_sen_phytomer",
                    "TT_del_phytomer",
                ],
                dtype=float,
            )

            self.phenT_tmp["id_phen"] = id_phen_list
            self.phenT_tmp["index_phytomer"] = index_phytomer_list

            phenT_tmp_grouped = self.phenT_tmp.groupby("id_phen")
            dynT_grouped = dynT_.groupby(["id_cohort", "N_phytomer_potential"])

            a_cohort_multiplicative_factor = 1 + params.A_COHORT_REDUCTION_FACTOR

            most_frequent_MS_a_cohort = dynT_["a_cohort"][dynT_.first_valid_index()]
            most_frequent_MS_TT_hs_0 = dynT_["TT_hs_0"][dynT_.first_valid_index()]
            most_frequent_MS_TT_start_MS_elongation = (
                decimal_elongated_internode_number / most_frequent_MS_a_cohort
            ) + most_frequent_MS_TT_hs_0

            for (id_cohort, N_phytomer_potential, id_phen), axeT_group in axeT_.groupby(
                ["id_cohort", "N_phytomer_potential", "id_phen"]
            ):
                phenT_tmp_group = phenT_tmp_grouped.get_group(id_phen)
                dynT_group = dynT_grouped.get_group((id_cohort, N_phytomer_potential))
                dynT_row = dynT_group.loc[
                    dynT_group[
                        dynT_group["id_axis"] == dynT_group["id_axis"].max()
                    ].first_valid_index()
                ]
                (
                    a_cohort_before_start_MS_elongation_1,
                    TT_hs_0,
                    TT_hs_break,
                    TT_flag_ligulation,
                ) = dynT_row[
                    ["a_cohort", "TT_hs_0", "TT_hs_break", "TT_flag_ligulation"]
                ]

                is_regressive = not bool(int(str(int(id_phen))[-1]))
                tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col = (
                    a_cohort_before_start_MS_elongation_1
                    * (most_frequent_MS_TT_start_MS_elongation - TT_hs_0)
                )
                tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip = (
                    tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col
                    - params.DELAIS_PHYLL_COL_TIP_NTH
                    / a_cohort_before_start_MS_elongation_1
                )

                if is_regressive:
                    a_cohort_after_start_MS_elongation_1 = (
                        a_cohort_before_start_MS_elongation_1
                        * a_cohort_multiplicative_factor
                    )
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col = (
                        tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col
                        - a_cohort_after_start_MS_elongation_1
                        * most_frequent_MS_TT_start_MS_elongation
                    )
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_tip = (
                        tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip
                        - a_cohort_after_start_MS_elongation_1
                        * most_frequent_MS_TT_start_MS_elongation
                    )  # TODO: useful ?
                else:
                    a_cohort_after_start_MS_elongation_1 = (
                        a_cohort_before_start_MS_elongation_1
                    )
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col = (
                        -a_cohort_before_start_MS_elongation_1 * TT_hs_0
                    )
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_tip = (
                        tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col
                        * params.DELAIS_PHYLL_COL_TIP_NTH
                        / a_cohort_before_start_MS_elongation_1
                    )  # TODO: useful ?

                if math.isnan(TT_hs_break):  # linear mode
                    HS_break = None
                    a_cohort_before_start_MS_elongation_2 = None
                    a_cohort_after_start_MS_elongation_2 = None
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col = None
                    tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_tip = None
                else:
                    HS_break = a_cohort_before_start_MS_elongation_1 * (
                        TT_hs_break - TT_hs_0
                    )
                    a_cohort_before_start_MS_elongation_2 = (
                        N_phytomer_potential - HS_break
                    ) / (TT_flag_ligulation - TT_hs_break)
                    if is_regressive:
                        a_cohort_after_start_MS_elongation_2 = (
                            a_cohort_before_start_MS_elongation_2
                            * a_cohort_multiplicative_factor
                        )
                        tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col = (
                            tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col
                            - a_cohort_after_start_MS_elongation_2
                            * most_frequent_MS_TT_start_MS_elongation
                        )
                        tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_tip = (
                            tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip
                            - a_cohort_after_start_MS_elongation_2
                            * most_frequent_MS_TT_start_MS_elongation
                        )  # TODO: useful ?
                    else:
                        a_cohort_after_start_MS_elongation_2 = (
                            a_cohort_before_start_MS_elongation_2
                        )
                        tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col = (
                            HS_break
                            - a_cohort_before_start_MS_elongation_2 * TT_hs_break
                        )
                        tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_tip = (
                            tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col
                            * params.DELAIS_PHYLL_COL_TIP_NTH
                            / a_cohort_before_start_MS_elongation_2
                        )  # TODO: useful ?

                # compute TT_col_phytomer
                self.phenT_tmp.loc[
                    phenT_tmp_group.index, "TT_col_phytomer"
                ] = phenT_tmp_group.loc[:, "TT_col_phytomer"].values[:] = (
                    phenT_tmp_group["index_phytomer"].apply(
                        _calculate_TT_col_phytomer,
                        args=(
                            HS_break,
                            TT_hs_0,
                            TT_hs_break,
                            tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col,
                            a_cohort_before_start_MS_elongation_1,
                            a_cohort_after_start_MS_elongation_1,
                            tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col,
                            a_cohort_before_start_MS_elongation_2,
                            a_cohort_after_start_MS_elongation_2,
                            tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col,
                            is_regressive,
                        ),
                    )
                )

                # compute TT_em_phytomer
                first_leaf_indexes = phenT_tmp_group.index[0:2]
                self.phenT_tmp.loc[
                    first_leaf_indexes, "TT_em_phytomer"
                ] = phenT_tmp_group.loc[first_leaf_indexes, "TT_em_phytomer"].values[
                    :
                ] = phenT_tmp_group.loc[
                    first_leaf_indexes, ["index_phytomer", "TT_col_phytomer"]
                ].apply(
                    _calculate_TT_em_phytomer,
                    axis=1,
                    args=(
                        TT_hs_break,
                        params.DELAIS_PHYLL_COL_TIP_1ST,
                        tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip,
                        a_cohort_before_start_MS_elongation_1,
                        a_cohort_after_start_MS_elongation_1,
                        a_cohort_before_start_MS_elongation_2,
                        a_cohort_after_start_MS_elongation_2,
                        is_regressive,
                    ),
                )

                try:
                    other_leaves_indexes = phenT_tmp_group.index.difference(
                        phenT_tmp_group.index[0:2]
                    )
                except AttributeError:  # backward compatibility with pandas < 0.16
                    other_leaves_indexes = (
                        phenT_tmp_group.index - phenT_tmp_group.index[0:2]
                    )
                if len(other_leaves_indexes) != 0:
                    self.phenT_tmp.loc[
                        other_leaves_indexes, "TT_em_phytomer"
                    ] = phenT_tmp_group.loc[
                        other_leaves_indexes, "TT_em_phytomer"
                    ].values[:] = phenT_tmp_group.loc[
                        other_leaves_indexes, ["index_phytomer", "TT_col_phytomer"]
                    ].apply(
                        _calculate_TT_em_phytomer,
                        axis=1,
                        args=(
                            TT_hs_break,
                            params.DELAIS_PHYLL_COL_TIP_NTH,
                            tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip,
                            a_cohort_before_start_MS_elongation_1,
                            a_cohort_after_start_MS_elongation_1,
                            a_cohort_before_start_MS_elongation_2,
                            a_cohort_after_start_MS_elongation_2,
                            is_regressive,
                        ),
                    )

                # compute TT_sen_phytomer
                id_axis, n0, n1, n2, t0, t1, a, c = dynT_row[
                    ["id_axis", "n0", "n1", "n2", "t0", "t1", "a", "c"]
                ]
                HS_1, HS_2, GL_2, GL_3, GL_4 = _calculate_HS_GL_polynomial(
                    HS_break,
                    id_axis,
                    a_cohort_before_start_MS_elongation_1,
                    TT_hs_0,
                    TT_hs_break,
                    TT_flag_ligulation,
                    n0,
                    n1,
                    n2,
                    t0,
                    t1,
                    a,
                    c,
                    a_cohort_before_start_MS_elongation_2,
                )

                self.phenT_tmp.loc[
                    phenT_tmp_group.index, "TT_sen_phytomer"
                ] = phenT_tmp_group.loc[:, "TT_sen_phytomer"].values[:] = (
                    phenT_tmp_group["index_phytomer"].apply(
                        _calculate_TT_sen_phytomer,
                        args=(
                            HS_break,
                            HS_1,
                            HS_2,
                            GL_2,
                            GL_3,
                            GL_4,
                            t0,
                            t1,
                            TT_flag_ligulation,
                            N_phytomer_potential,
                            is_regressive,
                        ),
                    )
                )

                # compute TT_del_phytomer
                self.phenT_tmp.loc[
                    phenT_tmp_group.index, "TT_del_phytomer"
                ] = phenT_tmp_group.loc[:, "TT_del_phytomer"].values[:] = (
                    _calculate_TT_del_phytomer(
                        a_cohort_before_start_MS_elongation_1,
                        phenT_tmp_group["TT_sen_phytomer"],
                    )
                )
        return self.phenT_tmp


_create_phenT_tmp = _CreatePhenTTmp()


def _create_phenT_abs(phenT_tmp, axeT_, dimT_):
    """
    Create the :ref:`phenT_abs <phenT_abs>` dataframe.
    """
    phenT_abs = phenT_tmp.copy()
    axeT_grouped = axeT_.groupby("id_phen")
    dimT_grouped = dimT_.groupby("id_dim")
    for id_phen, phenT_abs_group in phenT_abs.groupby("id_phen"):
        axeT_group = axeT_grouped.get_group(id_phen)
        id_dim = axeT_group["id_dim"][axeT_group.first_valid_index()]
        dimT_group = dimT_grouped.get_group(id_dim)
        non_zero_L_internode_group = dimT_group[dimT_group["L_internode"] != 0]
        if len(non_zero_L_internode_group) > 0:
            min_index_phytomer = non_zero_L_internode_group["index_phytomer"].min()
            indexes_to_ceil = phenT_abs_group[
                phenT_abs_group["index_phytomer"] >= min_index_phytomer
            ].index
            phenT_abs.loc[indexes_to_ceil, "TT_del_phytomer"] = params.TT_DEL_FHAUT
    return phenT_abs


def _calculate_TT_col_phytomer(
    index_phytomer,
    HS_break,
    TT_hs_0,
    TT_hs_break,
    tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col,
    a_cohort_before_start_MS_elongation_1,
    a_cohort_after_start_MS_elongation_1,
    tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col,
    a_cohort_before_start_MS_elongation_2,
    a_cohort_after_start_MS_elongation_2,
    tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col,
    is_regressive,
):
    if (
        math.isnan(TT_hs_break)
        or (index_phytomer + params.DELAIS_PHYLL_HS_COL_NTH) < HS_break
    ):  # 1st phase
        if (
            index_phytomer
            <= tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col
            or not is_regressive
        ):  # non regressive
            TT_col_phytomer = (
                index_phytomer + params.DELAIS_PHYLL_HS_COL_NTH
            ) / a_cohort_before_start_MS_elongation_1 + TT_hs_0
        else:  # regressive
            TT_col_phytomer = (
                index_phytomer
                + params.DELAIS_PHYLL_HS_COL_NTH
                - tiller_index_phytomer_y_intercept_after_start_MS_elongation_1_col
            ) / a_cohort_after_start_MS_elongation_1
    else:  # 2nd phase: bilinear mode only
        if (
            index_phytomer
            <= tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_col
            or not is_regressive
        ):  # non regressive
            TT_col_phytomer = (
                index_phytomer + params.DELAIS_PHYLL_HS_COL_NTH - HS_break
            ) / a_cohort_before_start_MS_elongation_2 + TT_hs_break
        else:  # regressive
            TT_col_phytomer = (
                index_phytomer
                + params.DELAIS_PHYLL_HS_COL_NTH
                - tiller_index_phytomer_y_intercept_after_start_MS_elongation_2_col
            ) / a_cohort_after_start_MS_elongation_2
    return TT_col_phytomer


def _calculate_TT_em_phytomer(
    phenT_subgroup,
    TT_hs_break,
    delais_phyll_col_tip,
    tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip,
    a_cohort_before_start_MS_elongation_1,
    a_cohort_after_start_MS_elongation_1,
    a_cohort_before_start_MS_elongation_2,
    a_cohort_after_start_MS_elongation_2,
    is_regressive,
):
    index_phytomer, TT_col_phytomer = (
        phenT_subgroup.index_phytomer,
        phenT_subgroup.TT_col_phytomer,
    )
    if math.isnan(TT_hs_break) or TT_col_phytomer < TT_hs_break:  # 1st phase
        if (
            index_phytomer
            <= tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip
            or not is_regressive
        ):  # non regressive
            TT_em_phytomer = (
                TT_col_phytomer
                - delais_phyll_col_tip / a_cohort_before_start_MS_elongation_1
            )
        else:  # regressive
            TT_em_phytomer = (
                TT_col_phytomer
                - delais_phyll_col_tip / a_cohort_after_start_MS_elongation_1
            )
    else:  # 2nd phase: bilinear mode only
        if (
            index_phytomer
            <= tiller_index_phytomer_at_most_frequent_MS_start_MS_elongation_tip
            or not is_regressive
        ):  # non regressive
            TT_em_phytomer = (
                TT_col_phytomer
                - delais_phyll_col_tip / a_cohort_before_start_MS_elongation_2
            )
            if TT_em_phytomer <= TT_hs_break:
                TT_em_phytomer = (
                    TT_hs_break
                    - (
                        delais_phyll_col_tip
                        - a_cohort_before_start_MS_elongation_2
                        * (TT_col_phytomer - TT_hs_break)
                    )
                    / a_cohort_before_start_MS_elongation_1
                )
        else:  # regressive
            TT_em_phytomer = (
                TT_col_phytomer
                - delais_phyll_col_tip / a_cohort_after_start_MS_elongation_2
            )
            if TT_em_phytomer <= TT_hs_break:
                TT_em_phytomer = (
                    TT_hs_break
                    - (
                        delais_phyll_col_tip
                        - a_cohort_after_start_MS_elongation_2
                        * (TT_col_phytomer - TT_hs_break)
                    )
                    / a_cohort_after_start_MS_elongation_1
                )
    return TT_em_phytomer


def _calculate_TT_sen_phytomer(
    index_phytomer,
    HS_break,
    HS_1,
    HS_2,
    GL_2,
    GL_3,
    GL_4,
    t0,
    t1,
    TT_flag_ligulation,
    N_phytomer_potential,
    is_regressive,
):
    # define HS according to index_phytomer
    if HS_break is None or index_phytomer < HS_break:  # linear mode
        HS = HS_1
    else:  # bilinear mode
        HS = HS_2
    GL_1 = HS

    if np.isnan(t0):  # 3 phases
        # Suppose TT_sen_phytomer is in the phase [0,t1].
        GL = GL_2
        SSI = HS - GL
        if is_regressive:
            SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
        if index_phytomer == 0:
            return t1
        else:
            # Find (SSI - index_phytomer) real root.
            SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
            if SSI_root_array.size != 0 and (
                SSI_root_array[0] <= t1 or np.allclose(SSI_root_array[0], t1)
            ):
                return SSI_root_array[0]

    else:  # 4 phases
        # Suppose TT_sen_phytomer is in the phase [0,t0].
        GL = GL_1
        SSI = HS - GL
        if is_regressive:
            SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
        if index_phytomer == 0:
            return t0
        else:
            # Find (SSI - index_phytomer) real root.
            SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
            if SSI_root_array.size != 0 and (
                SSI_root_array[0] <= t0 or np.allclose(SSI_root_array[0], t0)
            ):
                return SSI_root_array[0]
            else:
                # Suppose TT_sen_phytomer is in the phase ]t0,t1].
                GL = GL_2
                SSI = HS - GL
                if is_regressive:
                    SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
                # Find (SSI - index_phytomer) real root.
                SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
                if (
                    SSI_root_array.size != 0
                    and (SSI_root_array[0] > t0 or np.allclose(SSI_root_array[0], t0))
                    and (SSI_root_array[0] <= t1 or np.allclose(SSI_root_array[0], t1))
                ):
                    return SSI_root_array[0]

    # Suppose TT_sen_phytomer is in the phase ]t1,TT_flag_ligulation].
    GL = GL_3
    SSI = HS - GL
    if is_regressive:
        SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
    # Find (SSI - index_phytomer) real root.
    SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
    if (
        SSI_root_array.size != 0
        and (SSI_root_array[0] > t1 or np.allclose(SSI_root_array[0], t1))
        and (
            SSI_root_array[0] <= TT_flag_ligulation
            or np.allclose(SSI_root_array[0], TT_flag_ligulation)
        )
    ):
        return SSI_root_array[0]

    # TT_sen_phytomer must be in the phase ]TT_flag_ligulation,infinity[.
    GL = GL_4
    SSI = HS - GL
    if is_regressive:
        SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
    # Find (SSI - index_phytomer) real root.
    SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
    if SSI_root_array.size == 0 or (
        SSI_root_array[0] <= TT_flag_ligulation
        and not np.allclose(SSI_root_array[0], TT_flag_ligulation)
    ):
        raise Exception("ERROR !!!!! This shouldn't occurred")
    if HS(SSI_root_array[0]) > N_phytomer_potential:
        HS = np.poly1d([N_phytomer_potential])
    if GL(SSI_root_array[0]) < 0.0:
        GL = np.poly1d([0.0])
    SSI = HS - GL
    if is_regressive:
        SSI += params.MS_TO_REGRESSIVE_TILLERS_SENESCENCE_DELAY
    # Find (SSI - index_phytomer) real root again.
    SSI_root_array = tools.get_real_roots(SSI - index_phytomer)
    if SSI_root_array.size == 0 or (
        SSI_root_array[0] <= TT_flag_ligulation
        and not np.allclose(SSI_root_array[0], TT_flag_ligulation)
    ):
        raise Exception("ERROR !!!!! This shouldn't occurred")
    TT_sen_phytomer = SSI_root_array[0]

    return TT_sen_phytomer


def _calculate_HS_GL_polynomial(
    HS_break,
    id_axis,
    a_cohort_before_start_MS_elongation_1,
    TT_hs_0,
    TT_hs_break,
    TT_flag_ligulation,
    n0,
    n1,
    n2,
    t0,
    t1,
    a,
    c,
    a_cohort_before_start_MS_elongation_2,
):
    # define HS(TT)
    HS_1 = np.poly1d(
        [
            a_cohort_before_start_MS_elongation_1,
            -a_cohort_before_start_MS_elongation_1 * TT_hs_0,
        ]
    )  # index_phytomer < HS_break
    if HS_break is None:  # linear
        HS_2 = None  # index_phytomer >= HS_break
    else:  # bilinear
        HS_2 = np.poly1d(
            [
                a_cohort_before_start_MS_elongation_2,
                -a_cohort_before_start_MS_elongation_2 * TT_hs_break + HS_break,
            ]
        )  # index_phytomer >= HS_break

    # define GL(TT) for all phases except TT < t0 (because it depends on index_phytomer)
    if id_axis == "MS":
        GL_2 = np.poly1d([(n1 - n0) / (t1 - t0), n0 - t0 * (n1 - n0) / (t1 - t0)])
        GL_3 = np.poly1d(
            [
                (n2 - n1) / (TT_flag_ligulation - t1),
                n1 - t1 * (n2 - n1) / (TT_flag_ligulation - t1),
            ]
        )
    else:  # tillers
        if np.isnan(t0):  # only 3 phases
            GL_2 = np.poly1d([n1 / (t1 - TT_hs_0), n1 * TT_hs_0 / (TT_hs_0 - t1)])
        else:
            GL_2 = np.poly1d([(n1 - n0) / (t1 - t0), n1 - t1 * (n1 - n0) / (t1 - t0)])
        GL_3 = np.poly1d(
            [
                (n2 - n1) / (TT_flag_ligulation - t1),
                n2 - TT_flag_ligulation * (n2 - n1) / (TT_flag_ligulation - t1),
            ]
        )
    GL_4 = np.poly1d(
        [
            a,
            -3 * a * TT_flag_ligulation,
            3 * a * TT_flag_ligulation**2 + c,
            -a * TT_flag_ligulation**3 - c * TT_flag_ligulation + n2,
        ]
    )
    return HS_1, HS_2, GL_2, GL_3, GL_4


def _calculate_TT_del_phytomer(
    a_cohort_before_start_MS_elongation_1, TT_sen_phytomer_series
):
    TT_del_phytomer_series = (
        TT_sen_phytomer_series
        + params.DELAIS_PHYLL_SEN_DISP / a_cohort_before_start_MS_elongation_1
    )
    return TT_del_phytomer_series


class _CreatePhenTFirst:
    """
    Create the :ref:`phenT_first <phenT_first>` dataframe.
    """

    def __init__(self):
        self.phenT_first = None

    def __call__(self, phenT_tmp, force=True):
        if force or self.phenT_first is None:
            if not (
                phenT_tmp.count().max()
                == phenT_tmp.count().min()
                == phenT_tmp.index.size
            ):
                raise tools.InputError("phenT_tmp contains NA values")
            self.phenT_first = phenT_tmp.loc[
                phenT_tmp.index.map(
                    lambda idx: True if phenT_tmp["index_phytomer"][idx] == 1 else False
                )
            ]
        return self.phenT_first


_create_phenT_first = _CreatePhenTFirst()


def _create_phenT(phenT_abs, phenT_first):
    """
    Create the :ref:`phenT <phenT>` dataframe.
    """
    if not all(phenT_abs.notnull()):
        raise tools.InputError("phenT_abs contains NA values")

    # define TT_*_phytomer_1 and merge in tmp
    tmp_first = phenT_first.drop("index_phytomer", axis=1)
    stades = ("em", "col", "sen", "del")
    tmp_first = tmp_first.rename(
        columns={
            "_".join(("TT", k, "phytomer")): "_".join(("TT", k, "phytomer", "1"))
            for k in stades
        }
    )
    tmp = pd.merge(phenT_abs, tmp_first, on="id_phen")

    # build phenT_
    phenT_ = pd.DataFrame(
        index=phenT_abs.index,
        columns=[
            "id_phen",
            "index_phytomer",
            "dTT_em_phytomer",
            "dTT_col_phytomer",
            "dTT_sen_phytomer",
            "dTT_del_phytomer",
        ],
        dtype=float,
    )
    phenT_["id_phen"] = phenT_abs["id_phen"]
    phenT_["index_phytomer"] = tmp["index_phytomer"]
    for w in stades:
        try:
            phenT_["_".join(("dTT", w, "phytomer"))] = (
                tmp["_".join(("TT", w, "phytomer"))]
                - tmp["_".join(("TT", w, "phytomer", "1"))]
            ).tolist()
        except ValueError:
            pass
    return phenT_


def _create_HS_GL_SSI_T(axeT_, dynT_):
    """
    Create the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>` dataframe.
    """
    HS_GL_SSI_dynamic = pd.DataFrame(columns=["id_phen", "TT", "HS", "GL", "SSI"])

    dynT_grouped = dynT_.groupby(["id_cohort", "N_phytomer_potential"])
    for (id_cohort, N_phytomer_potential, id_phen), axeT_group in axeT_.groupby(
        ["id_cohort", "N_phytomer_potential", "id_phen"]
    ):
        dynT_group = dynT_grouped.get_group((id_cohort, N_phytomer_potential))
        dynT_row = dynT_group.loc[dynT_group["cardinality"].idxmax()]
        (
            id_axis,
            a_cohort_before_start_MS_elongation_1,
            TT_hs_0,
            TT_hs_break,
            TT_flag_ligulation,
            n0,
            n1,
            n2,
            t0,
            t1,
            a,
            c,
        ) = dynT_row[
            [
                "id_axis",
                "a_cohort",
                "TT_hs_0",
                "TT_hs_break",
                "TT_flag_ligulation",
                "n0",
                "n1",
                "n2",
                "t0",
                "t1",
                "a",
                "c",
            ]
        ]

        if math.isnan(TT_hs_break):  # linear mode
            HS_break = None
            a_cohort_before_start_MS_elongation_2 = None
        else:  # bilinear mode
            HS_break = a_cohort_before_start_MS_elongation_1 * (TT_hs_break - TT_hs_0)
            a_cohort_before_start_MS_elongation_2 = (
                N_phytomer_potential - HS_break
            ) / (TT_flag_ligulation - TT_hs_break)

        HS_1, HS_2, GL_2, GL_3, GL_4 = _calculate_HS_GL_polynomial(
            HS_break,
            id_axis,
            a_cohort_before_start_MS_elongation_1,
            TT_hs_0,
            TT_hs_break,
            TT_flag_ligulation,
            n0,
            n1,
            n2,
            t0,
            t1,
            a,
            c,
            a_cohort_before_start_MS_elongation_2,
        )

        t1 = int(round(t1))
        TT_flag_ligulation = int(round(TT_flag_ligulation))

        if np.isnan(t0):  # 3 phases
            TT_2 = np.arange(0, t1)
            TT_3 = np.arange(t1, TT_flag_ligulation)
            TT_4 = np.arange(TT_flag_ligulation, params.TT_DEL_FHAUT)

            if math.isnan(TT_hs_break):  # linear mode
                HS_1_TT_2 = np.clip(HS_1(TT_2), 0.0, N_phytomer_potential)
                HS_1_TT_3 = np.clip(HS_1(TT_3), 0.0, N_phytomer_potential)
                HS_1_TT_4 = np.clip(HS_1(TT_4), 0.0, N_phytomer_potential)

                GL_2_TT_2 = np.clip(GL_2(TT_2), 0.0, 1000.0)
                GL_3_TT_3 = np.clip(GL_3(TT_3), 0.0, 1000.0)
                GL_4_TT_4 = np.clip(GL_4(TT_4), 0.0, 1000.0)

                SSI_2_TT_2 = HS_1_TT_2 - GL_2_TT_2
                SSI_3_TT_3 = HS_1_TT_3 - GL_3_TT_3
                SSI_4_TT_4 = HS_1_TT_4 - GL_4_TT_4

                HS_GL_SSI_dynamic_group = pd.DataFrame(
                    index=np.arange(params.TT_DEL_FHAUT),
                    columns=HS_GL_SSI_dynamic.columns,
                )
                HS_GL_SSI_dynamic_group["id_phen"] = np.repeat(
                    id_phen, HS_GL_SSI_dynamic_group.index.size
                )
                HS_GL_SSI_dynamic_group["TT"] = pd.Series(
                    np.concatenate((TT_2, TT_3, TT_4))
                )
                HS_GL_SSI_dynamic_group["HS"] = pd.Series(
                    np.concatenate((HS_1_TT_2, HS_1_TT_3, HS_1_TT_4))
                )
                HS_GL_SSI_dynamic_group["GL"] = pd.Series(
                    np.concatenate((GL_2_TT_2, GL_3_TT_3, GL_4_TT_4))
                )
                HS_GL_SSI_dynamic_group["SSI"] = pd.Series(
                    np.concatenate((SSI_2_TT_2, SSI_3_TT_3, SSI_4_TT_4))
                )

            else:  # bilinear mode
                TT_hs_break = int(round(TT_hs_break))
                TT_2_1, TT_3_1, TT_4_1 = TT_2, TT_3, TT_4
                TT_2_2 = TT_3_2 = TT_4_2 = []
                if TT_hs_break <= t1:
                    TT_2_1 = np.arange(0, TT_hs_break)
                    TT_2_2 = np.arange(TT_hs_break, t1)
                elif TT_hs_break <= TT_flag_ligulation:
                    TT_3_1 = np.arange(t1, TT_hs_break)
                    TT_3_2 = np.arange(TT_hs_break, TT_flag_ligulation)
                else:
                    TT_4_1 = np.arange(TT_flag_ligulation, TT_hs_break)
                    TT_4_2 = np.arange(TT_hs_break, params.TT_DEL_FHAUT)

                HS_1_TT_2_1 = np.clip(HS_1(TT_2_1), 0.0, N_phytomer_potential)
                HS_1_TT_3_1 = np.clip(HS_1(TT_3_1), 0.0, N_phytomer_potential)
                HS_1_TT_4_1 = np.clip(HS_1(TT_4_1), 0.0, N_phytomer_potential)
                HS_1_TT_2_2 = np.clip(HS_1(TT_2_2), 0.0, N_phytomer_potential)
                HS_1_TT_3_2 = np.clip(HS_1(TT_3_2), 0.0, N_phytomer_potential)
                HS_1_TT_4_2 = np.clip(HS_1(TT_4_2), 0.0, N_phytomer_potential)

                GL_2_TT_2_1 = np.clip(GL_2(TT_2_1), 0.0, 1000.0)
                GL_3_TT_3_1 = np.clip(GL_3(TT_3_1), 0.0, 1000.0)
                GL_4_TT_4_1 = np.clip(GL_4(TT_4_1), 0.0, 1000.0)
                GL_2_TT_2_2 = np.clip(GL_2(TT_2_2), 0.0, 1000.0)
                GL_3_TT_3_2 = np.clip(GL_3(TT_3_2), 0.0, 1000.0)
                GL_4_TT_4_2 = np.clip(GL_4(TT_4_2), 0.0, 1000.0)

                SSI_2_TT_2_1 = HS_1_TT_2_1 - GL_2_TT_2_1
                SSI_3_TT_3_1 = HS_1_TT_3_1 - GL_3_TT_3_1
                SSI_4_TT_4_1 = HS_1_TT_4_1 - GL_4_TT_4_1
                SSI_2_TT_2_2 = HS_1_TT_2_2 - GL_2_TT_2_2
                SSI_3_TT_3_2 = HS_1_TT_3_2 - GL_3_TT_3_2
                SSI_4_TT_4_2 = HS_1_TT_4_2 - GL_4_TT_4_2

                HS_GL_SSI_dynamic_group = pd.DataFrame(
                    index=np.arange(params.TT_DEL_FHAUT),
                    columns=HS_GL_SSI_dynamic.columns,
                )
                HS_GL_SSI_dynamic_group["id_phen"] = np.repeat(
                    id_phen, HS_GL_SSI_dynamic_group.index.size
                )
                HS_GL_SSI_dynamic_group["TT"] = pd.Series(
                    np.concatenate((TT_2_1, TT_2_2, TT_3_1, TT_3_2, TT_4_1, TT_4_2))
                )
                HS_GL_SSI_dynamic_group["HS"] = pd.Series(
                    np.concatenate(
                        (
                            HS_1_TT_2_1,
                            HS_1_TT_2_2,
                            HS_1_TT_3_1,
                            HS_1_TT_3_2,
                            HS_1_TT_4_1,
                            HS_1_TT_4_2,
                        )
                    )
                )
                HS_GL_SSI_dynamic_group["GL"] = pd.Series(
                    np.concatenate(
                        (
                            GL_2_TT_2_1,
                            GL_2_TT_2_2,
                            GL_3_TT_3_1,
                            GL_3_TT_3_2,
                            GL_4_TT_4_1,
                            GL_4_TT_4_2,
                        )
                    )
                )
                HS_GL_SSI_dynamic_group["SSI"] = pd.Series(
                    np.concatenate(
                        (
                            SSI_2_TT_2_1,
                            SSI_2_TT_2_2,
                            SSI_3_TT_3_1,
                            SSI_3_TT_3_2,
                            SSI_4_TT_4_1,
                            SSI_4_TT_4_2,
                        )
                    )
                )

            HS_GL_SSI_dynamic = pd.concat(
                [HS_GL_SSI_dynamic, HS_GL_SSI_dynamic_group], ignore_index=True
            )

        else:  # 4 phases
            t0 = int(round(t0))

            TT_1 = np.arange(0, t0)
            TT_2 = np.arange(t0, t1)
            TT_3 = np.arange(t1, TT_flag_ligulation)
            TT_4 = np.arange(TT_flag_ligulation, params.TT_DEL_FHAUT)

            if math.isnan(TT_hs_break):  # linear mode
                HS_1_TT_1 = np.clip(HS_1(TT_1), 0.0, N_phytomer_potential)
                HS_1_TT_2 = np.clip(HS_1(TT_2), 0.0, N_phytomer_potential)
                HS_1_TT_3 = np.clip(HS_1(TT_3), 0.0, N_phytomer_potential)
                HS_1_TT_4 = np.clip(HS_1(TT_4), 0.0, N_phytomer_potential)

                GL_1_TT_1 = np.clip(HS_1(TT_1), 0.0, 1000.0)
                GL_2_TT_2 = np.clip(GL_2(TT_2), 0.0, 1000.0)
                GL_3_TT_3 = np.clip(GL_3(TT_3), 0.0, 1000.0)
                GL_4_TT_4 = np.clip(GL_4(TT_4), 0.0, 1000.0)

                SSI_1_TT_1 = HS_1_TT_1 - GL_1_TT_1
                SSI_2_TT_2 = HS_1_TT_2 - GL_2_TT_2
                SSI_3_TT_3 = HS_1_TT_3 - GL_3_TT_3
                SSI_4_TT_4 = HS_1_TT_4 - GL_4_TT_4

                HS_GL_SSI_dynamic_group = pd.DataFrame(
                    index=np.arange(params.TT_DEL_FHAUT),
                    columns=HS_GL_SSI_dynamic.columns,
                )
                HS_GL_SSI_dynamic_group["id_phen"] = np.repeat(
                    id_phen, HS_GL_SSI_dynamic_group.index.size
                )
                HS_GL_SSI_dynamic_group["TT"] = pd.Series(
                    np.concatenate((TT_1, TT_2, TT_3, TT_4))
                )
                HS_GL_SSI_dynamic_group["HS"] = pd.Series(
                    np.concatenate((HS_1_TT_1, HS_1_TT_2, HS_1_TT_3, HS_1_TT_4))
                )
                HS_GL_SSI_dynamic_group["GL"] = pd.Series(
                    np.concatenate((GL_1_TT_1, GL_2_TT_2, GL_3_TT_3, GL_4_TT_4))
                )
                HS_GL_SSI_dynamic_group["SSI"] = pd.Series(
                    np.concatenate((SSI_1_TT_1, SSI_2_TT_2, SSI_3_TT_3, SSI_4_TT_4))
                )

            else:  # bilinear mode
                TT_hs_break = int(round(TT_hs_break))
                TT_1_1, TT_2_1, TT_3_1, TT_4_1 = TT_1, TT_2, TT_3, TT_4
                TT_1_2 = TT_2_2 = TT_3_2 = TT_4_2 = []
                if TT_hs_break <= t0:
                    TT_1_1 = np.arange(0, TT_hs_break)
                    TT_1_2 = np.arange(TT_hs_break, t0)
                elif TT_hs_break <= t1:
                    TT_2_1 = np.arange(t0, TT_hs_break)
                    TT_2_2 = np.arange(TT_hs_break, t1)
                elif TT_hs_break <= TT_flag_ligulation:
                    TT_3_1 = np.arange(t1, TT_hs_break)
                    TT_3_2 = np.arange(TT_hs_break, TT_flag_ligulation)
                else:
                    TT_4_1 = np.arange(TT_flag_ligulation, TT_hs_break)
                    TT_4_2 = np.arange(TT_hs_break, params.TT_DEL_FHAUT)

                HS_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, N_phytomer_potential)
                HS_1_TT_2_1 = np.clip(HS_1(TT_2_1), 0.0, N_phytomer_potential)
                HS_1_TT_3_1 = np.clip(HS_1(TT_3_1), 0.0, N_phytomer_potential)
                HS_1_TT_4_1 = np.clip(HS_1(TT_4_1), 0.0, N_phytomer_potential)
                HS_1_TT_1_2 = np.clip(HS_1(TT_1_2), 0.0, N_phytomer_potential)
                HS_1_TT_2_2 = np.clip(HS_1(TT_2_2), 0.0, N_phytomer_potential)
                HS_1_TT_3_2 = np.clip(HS_1(TT_3_2), 0.0, N_phytomer_potential)
                HS_1_TT_4_2 = np.clip(HS_1(TT_4_2), 0.0, N_phytomer_potential)

                GL_1_TT_1_1 = np.clip(HS_1(TT_1_1), 0.0, 1000.0)
                GL_2_TT_2_1 = np.clip(GL_2(TT_2_1), 0.0, 1000.0)
                GL_3_TT_3_1 = np.clip(GL_3(TT_3_1), 0.0, 1000.0)
                GL_4_TT_4_1 = np.clip(GL_4(TT_4_1), 0.0, 1000.0)
                GL_1_TT_1_2 = np.clip(HS_2(TT_1_2), 0.0, 1000.0)
                GL_2_TT_2_2 = np.clip(GL_2(TT_2_2), 0.0, 1000.0)
                GL_3_TT_3_2 = np.clip(GL_3(TT_3_2), 0.0, 1000.0)
                GL_4_TT_4_2 = np.clip(GL_4(TT_4_2), 0.0, 1000.0)

                SSI_1_TT_1_1 = HS_1_TT_1_1 - GL_1_TT_1_1
                SSI_2_TT_2_1 = HS_1_TT_2_1 - GL_2_TT_2_1
                SSI_3_TT_3_1 = HS_1_TT_3_1 - GL_3_TT_3_1
                SSI_4_TT_4_1 = HS_1_TT_4_1 - GL_4_TT_4_1
                SSI_1_TT_1_2 = HS_1_TT_1_2 - GL_1_TT_1_2
                SSI_2_TT_2_2 = HS_1_TT_2_2 - GL_2_TT_2_2
                SSI_3_TT_3_2 = HS_1_TT_3_2 - GL_3_TT_3_2
                SSI_4_TT_4_2 = HS_1_TT_4_2 - GL_4_TT_4_2

                HS_GL_SSI_dynamic_group = pd.DataFrame(
                    index=np.arange(params.TT_DEL_FHAUT),
                    columns=HS_GL_SSI_dynamic.columns,
                )
                HS_GL_SSI_dynamic_group["id_phen"] = np.repeat(
                    id_phen, HS_GL_SSI_dynamic_group.index.size
                )
                HS_GL_SSI_dynamic_group["TT"] = pd.Series(
                    np.concatenate(
                        (TT_1_1, TT_1_2, TT_2_1, TT_2_2, TT_3_1, TT_3_2, TT_4_1, TT_4_2)
                    )
                )
                HS_GL_SSI_dynamic_group["HS"] = pd.Series(
                    np.concatenate(
                        (
                            HS_1_TT_1_1,
                            HS_1_TT_1_2,
                            HS_1_TT_2_1,
                            HS_1_TT_2_2,
                            HS_1_TT_3_1,
                            HS_1_TT_3_2,
                            HS_1_TT_4_1,
                            HS_1_TT_4_2,
                        )
                    )
                )
                HS_GL_SSI_dynamic_group["GL"] = pd.Series(
                    np.concatenate(
                        (
                            GL_1_TT_1_1,
                            GL_1_TT_1_2,
                            GL_2_TT_2_1,
                            GL_2_TT_2_2,
                            GL_3_TT_3_1,
                            GL_3_TT_3_2,
                            GL_4_TT_4_1,
                            GL_4_TT_4_2,
                        )
                    )
                )
                HS_GL_SSI_dynamic_group["SSI"] = pd.Series(
                    np.concatenate(
                        (
                            SSI_1_TT_1_1,
                            SSI_1_TT_1_2,
                            SSI_2_TT_2_1,
                            SSI_2_TT_2_2,
                            SSI_3_TT_3_1,
                            SSI_3_TT_3_2,
                            SSI_4_TT_4_1,
                            SSI_4_TT_4_2,
                        )
                    )
                )

            HS_GL_SSI_dynamic = pd.concat(
                [HS_GL_SSI_dynamic if not HS_GL_SSI_dynamic.empty else None, HS_GL_SSI_dynamic_group], ignore_index=True
            )

    HS_GL_SSI_dynamic[["id_phen", "TT"]] = HS_GL_SSI_dynamic[["id_phen", "TT"]].astype(
        int
    )

    return HS_GL_SSI_dynamic


class _PrepareMergeWithUserTables:
    """
    Prepare the merge of:
        * :ref:`dynT_user` and dynT_tmp,
        * :ref:`dimT_user` and dimT_tmp.
    """

    def __init__(self):
        self.most_frequent_dynT_tmp = None
        self.most_frequent_dynT_tmp_grouped = None

    def __call__(self, dynT_tmp, force=True):
        if (
            force
            or self.most_frequent_dynT_tmp is None
            or self.most_frequent_dynT_tmp_grouped is None
        ):
            self.most_frequent_dynT_tmp = pd.DataFrame(columns=dynT_tmp.columns)
            for id_axis, dynT_tmp_group in dynT_tmp.groupby("id_axis"):
                idxmax = dynT_tmp_group["cardinality"].idxmax()
                self.most_frequent_dynT_tmp = pd.concat(
                    [self.most_frequent_dynT_tmp, dynT_tmp_group.loc[idxmax:idxmax]],
                    ignore_index=True,
                )
            self.most_frequent_dynT_tmp_grouped = self.most_frequent_dynT_tmp.groupby(
                ["id_axis", "N_phytomer_potential"]
            )
        return self.most_frequent_dynT_tmp, self.most_frequent_dynT_tmp_grouped


_prepare_merge_with_user_tables = _PrepareMergeWithUserTables()


def _merge_dynT_tmp_and_dynT_user(
    dynT_tmp, dynT_user, dynT_user_completeness, TT_hs_break
):
    """
    Merge :ref:`dynT_user` and dynT_tmp.
    """
    most_frequent_dynT_tmp, most_frequent_dynT_tmp_grouped = (
        _prepare_merge_with_user_tables(dynT_tmp)
    )
    dynT_tmp_merged = pd.DataFrame(dynT_tmp)
    if dynT_user_completeness == DataCompleteness.MIN:
        dynT_user_grouped = dynT_user.groupby("id_axis")
        MS_dynT_user = dynT_user_grouped.get_group("MS")
        MS_TT_flag_ligulation = MS_dynT_user["TT_flag_ligulation"][
            MS_dynT_user.first_valid_index()
        ]
        for id_axis, dynT_tmp_group in dynT_tmp_merged.groupby("id_axis"):
            if not dynT_user["id_axis"].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dynT_user_group = dynT_user_grouped.get_group(id_axis)
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            if id_axis == "MS":
                for column in [
                    "a_cohort",
                    "TT_hs_0",
                    "TT_flag_ligulation",
                    "n0",
                    "n1",
                    "n2",
                ]:
                    dynT_tmp_merged.loc[index_to_set, column] = dynT_user[column][
                        index_to_get
                    ]
            dynT_tmp_merged.loc[index_to_set, "TT_hs_break"] = TT_hs_break
            current_TT_flag_ligulation = dynT_user_group["TT_flag_ligulation"][
                index_to_get
            ]
            current_dTT_MS_cohort = current_TT_flag_ligulation - MS_TT_flag_ligulation
            dynT_tmp_merged.loc[index_to_set, "dTT_MS_cohort"] = current_dTT_MS_cohort
    elif dynT_user_completeness == DataCompleteness.SHORT:
        dynT_user_grouped = dynT_user.groupby("id_axis")
        MS_dynT_user = dynT_user_grouped.get_group("MS")
        MS_TT_flag_ligulation = MS_dynT_user["TT_flag_ligulation"][
            MS_dynT_user.first_valid_index()
        ]
        for id_axis, dynT_tmp_group in dynT_tmp_merged.groupby("id_axis"):
            if not dynT_user["id_axis"].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dynT_user_group = dynT_user_grouped.get_group(id_axis)
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            for column in [
                "a_cohort",
                "TT_hs_0",
                "TT_flag_ligulation",
                "n0",
                "n1",
                "n2",
            ]:
                dynT_tmp_merged.loc[index_to_set, column] = dynT_user[column][
                    index_to_get
                ]
            dynT_tmp_merged.loc[index_to_set, "TT_hs_break"] = TT_hs_break
            current_TT_flag_ligulation = dynT_user_group["TT_flag_ligulation"][
                index_to_get
            ]
            current_dTT_MS_cohort = current_TT_flag_ligulation - MS_TT_flag_ligulation
            dynT_tmp_merged.loc[index_to_set, "dTT_MS_cohort"] = current_dTT_MS_cohort
    elif dynT_user_completeness == DataCompleteness.FULL:
        dynT_user_grouped = dynT_user.groupby(["id_axis", "N_phytomer_potential"])
        MS_N_phytomer_potential = dynT_tmp_merged["N_phytomer_potential"][
            dynT_tmp_merged[dynT_tmp_merged["id_axis"] == "MS"]["cardinality"].idxmax()
        ]
        MS_dynT_user = dynT_user_grouped.get_group(("MS", MS_N_phytomer_potential))
        MS_TT_flag_ligulation = MS_dynT_user["TT_flag_ligulation"][
            MS_dynT_user.first_valid_index()
        ]
        for (id_axis, N_phytomer_potential), dynT_tmp_group in dynT_tmp_merged.groupby(
            ["id_axis", "N_phytomer_potential"]
        ):
            if not dynT_user["id_axis"].isin([id_axis]).any():
                if (
                    id_axis,
                    N_phytomer_potential,
                ) in most_frequent_dynT_tmp_grouped.groups:
                    id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
                else:
                    continue
            if (id_axis, N_phytomer_potential) not in dynT_user_grouped.groups:
                if (
                    id_axis,
                    N_phytomer_potential,
                ) in most_frequent_dynT_tmp_grouped.groups:
                    raise tools.InputError(
                        "Dynamic of %s not documented"
                        % ((id_axis, N_phytomer_potential),)
                    )
                else:
                    most_frequent_dynT_tmp_id_axis = most_frequent_dynT_tmp[
                        most_frequent_dynT_tmp["id_axis"] == id_axis
                    ]
                    N_phytomer_potential = most_frequent_dynT_tmp_id_axis[
                        "N_phytomer_potential"
                    ][most_frequent_dynT_tmp_id_axis.first_valid_index()]
                    if (id_axis, N_phytomer_potential) not in dynT_user_grouped.groups:
                        raise tools.InputError(
                            "Dynamic of %s not documented"
                            % ((id_axis, N_phytomer_potential),)
                        )
            dynT_user_group = dynT_user_grouped.get_group(
                (id_axis, N_phytomer_potential)
            )
            index_to_get = dynT_user_group.index[0]
            index_to_set = dynT_tmp_group.index[0]
            for column in [
                "a_cohort",
                "TT_hs_0",
                "TT_flag_ligulation",
                "n0",
                "n1",
                "n2",
            ]:
                dynT_tmp_merged.loc[index_to_set, column] = dynT_user[column][
                    index_to_get
                ]
            dynT_tmp_merged.loc[index_to_set, "TT_hs_break"] = TT_hs_break
            current_TT_flag_ligulation = dynT_user_group["TT_flag_ligulation"][
                index_to_get
            ]
            current_dTT_MS_cohort = current_TT_flag_ligulation - MS_TT_flag_ligulation
            dynT_tmp_merged.loc[index_to_set, "dTT_MS_cohort"] = current_dTT_MS_cohort

    return dynT_tmp_merged


def _merge_dimT_tmp_and_dimT_user(
    dynT_tmp_merged, dimT_user, dimT_user_completeness, dimT_tmp
):
    """
    Merge :ref:`dimT_user` and dimT_tmp.
    """
    most_frequent_dynT_tmp, most_frequent_dynT_tmp_grouped = (
        _prepare_merge_with_user_tables(dynT_tmp_merged, force=False)
    )
    dimT_tmp_merged = pd.DataFrame(dimT_tmp)
    organ_dim_list = [
        "L_blade",
        "W_blade",
        "L_sheath",
        "W_sheath",
        "L_internode",
        "W_internode",
    ]
    if dimT_user_completeness == DataCompleteness.MIN:
        MS_dynT_tmp = dynT_tmp_merged[dynT_tmp_merged["id_axis"] == "MS"]
        MS_most_frequent_N_phytomer_potential = MS_dynT_tmp["N_phytomer_potential"][
            MS_dynT_tmp["cardinality"].idxmax()
        ]
        dimT_tmp_grouped = dimT_tmp_merged.groupby(["id_axis", "N_phytomer_potential"])
        dimT_tmp_indexes_to_set = dimT_tmp_grouped.groups[
            ("MS", MS_most_frequent_N_phytomer_potential)
        ]
        N_phytomer_potential_to_set = len(dimT_tmp_indexes_to_set)
        max_available_MS_N_phytomer_potential = dimT_user["index_phytomer"].max()
        if N_phytomer_potential_to_set > max_available_MS_N_phytomer_potential:
            raise tools.InputError(
                "Dimensions of index_phytomer=%s not documented"
                % N_phytomer_potential_to_set
            )
        dimT_user_indexes_to_get = list(range(len(dimT_tmp_indexes_to_set)))
        for organ_dim in organ_dim_list:
            dimT_tmp_merged.loc[dimT_tmp_indexes_to_set, organ_dim] = dimT_user[
                organ_dim
            ][dimT_user_indexes_to_get].values
    elif dimT_user_completeness == DataCompleteness.SHORT:
        dimT_user_grouped = dimT_user.groupby("id_axis")
        dynT_tmp_grouped = dynT_tmp_merged.groupby("id_axis")
        for id_axis, dimT_tmp_group in dimT_tmp_merged.groupby("id_axis"):
            dynT_tmp_group = dynT_tmp_grouped.get_group(id_axis)
            N_phytomer_potential = dynT_tmp_group["N_phytomer_potential"][
                dynT_tmp_group["cardinality"].idxmax()
            ]
            indexes_to_set = dimT_tmp_group[
                dimT_tmp_group["N_phytomer_potential"] == N_phytomer_potential
            ].index
            if not dimT_user["id_axis"].isin([id_axis]).any():
                id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
            dimT_user_group = dimT_user_grouped.get_group(id_axis)
            max_available_N_phytomer_potential = dimT_user_group["index_phytomer"].max()
            N_phytomer_potential_to_set = len(indexes_to_set)
            if N_phytomer_potential_to_set > max_available_N_phytomer_potential:
                raise tools.InputError(
                    "Dimensions of %s not documented"
                    % ((id_axis, N_phytomer_potential_to_set),)
                )
            indexes_to_get = dimT_user_group.index[:N_phytomer_potential_to_set]
            for organ_dim in organ_dim_list:
                dimT_tmp_merged.loc[indexes_to_set, organ_dim] = dimT_user[organ_dim][
                    indexes_to_get
                ].values
    elif dimT_user_completeness == DataCompleteness.FULL:
        dimT_user_grouped = dimT_user.groupby(["id_axis", "N_phytomer_potential"])
        for (id_axis, N_phytomer_potential), dimT_tmp_group in dimT_tmp_merged.groupby(
            ["id_axis", "N_phytomer_potential"]
        ):
            if not dimT_user["id_axis"].isin([id_axis]).any():
                if (
                    id_axis,
                    N_phytomer_potential,
                ) in most_frequent_dynT_tmp_grouped.groups:
                    id_axis = tools.get_primary_axis(id_axis, params.FIRST_CHILD_DELAY)
                else:
                    continue
            indexes_to_set = dimT_tmp_group.index
            if (id_axis, N_phytomer_potential) not in dimT_user_grouped.groups:
                raise tools.InputError(
                    "Dimensions of %s not documented"
                    % ((id_axis, N_phytomer_potential),)
                )
            indexes_to_get = dimT_user_grouped.get_group(
                (id_axis, N_phytomer_potential)
            ).index
            for organ_dim in organ_dim_list:
                dimT_tmp_merged.loc[indexes_to_set, organ_dim] = dimT_user[organ_dim][
                    indexes_to_get
                ].values

    return dimT_tmp_merged
