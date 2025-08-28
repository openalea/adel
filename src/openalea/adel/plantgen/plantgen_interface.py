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
Front end for the generation of the input data expected by ADEL. User should
look at this module first. One can then look at the other modules of :mod:`openalea.adel.plantgen`
for additional information.

Authors: M. Abichou, B. Andrieu, C. Chambon
"""

import numpy as np
import pandas as pd
import warnings

from openalea.adel.plantgen import plantgen_core, tools, params

warnings.simplefilter("always", tools.InputWarning)


def gen_adel_input_data(
    dynT_user,
    dimT_user,
    plants_number=100,
    plants_density=250,
    decide_child_axis_probabilities={
        "T0": 0.0,
        "T1": 0.900,
        "T2": 0.983,
        "T3": 0.817,
        "T4": 0.117,
    },
    MS_leaves_number_probabilities={
        "10": 0.145,
        "11": 0.818,
        "12": 0.037,
        "13": 0.0,
        "14": 0.0,
    },
    ears_density=500,
    GL_number={1117.0: 5.6, 1212.1: 5.4, 1368.7: 4.9, 1686.8: 2.4, 1880.0: 0.0},
    delais_TT_stop_del_axis=600,
    TT_hs_break=None,
    inner_params={},
    axeT_user=None,
    TT_regression_start_user=None,
    TT_t1_user=None,
):
    """
    Create the dataframes which contain the plant data to be used as input for
    generating plot with ADEL, and some other dataframes for debugging purpose.
    Also create a dictionary which stores the values of the arguments of
    :func:`gen_adel_input_data <openalea.adel.plantgen.plantgen_interface.gen_adel_input_data>`.
    This dictionary is aimed to log the configuration used for the construction.

    See :ref:`adel_input` for a description of the input tables expected by ADEL,
    and :ref:`plantgen` for a description of the dataframes created for debug.

    Different degrees of completeness of data provided by the user are acceptable.
    See :ref:`levels_of_completeness` for more details.

    The dataframes are created as follows:
        * initialize the axes, calling :func:`init_axes <openalea.adel.plantgen.plantgen_core.init_axes>`,
        * define phenology functions, calling :func:`phenology_functions <openalea.adel.plantgen.plantgen_core.phenology_functions>`,
        * construct plants structure, calling :func:`plants_structure <openalea.adel.plantgen.plantgen_core.plants_structure>`,
        * calculate organs dimensions, calling :func:`organs_dimensions <openalea.adel.plantgen.plantgen_core.organs_dimensions>`,
        * calculate the phenology of the axes, calling :func:`axes_phenology <openalea.adel.plantgen.plantgen_core.axes_phenology>`.

    These tables are returned to be used as ADEL input:
        * the :ref:`axeT <axeT>`,
        * the :ref:`dimT <dimT>`,
        * the :ref:`phenT <phenT>`.

    These tables are returned for debugging purpose:
        * the :ref:`tilleringT <tilleringT>`,
        * the :ref:`cardinalityT <cardinalityT>`,
        * the :ref:`phenT_abs <phenT_abs>`,
        * the :ref:`dynT <dynT>`,
        * the :ref:`phenT_first <phenT_first>`,
        * the :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`,

    :Parameters:

        - `dynT_user` (:class:`pandas.DataFrame`) - the leaf dynamic
          parameters set by the user.

        - `dimT_user` (:class:`pandas.DataFrame`) - the dimensions of the organs
          set by the user.

        - `plants_number` (:class:`int`) - the number of plants to be generated.

        - `plants_density` (:class:`int`) - the number of plants that are present
          after loss due to bad emergence, early death..., per square meter.

        - `decide_child_axis_probabilities` (:class:`dict` of :class:`str`::class:`float`) -
          for each child cohort the probability of emergence of an axis when the parent
          axis is present. The keys are the identifiers of the child axes ('T0', 'T1',
          'T2', ...) and the values are the probabilities.

        - `MS_leaves_number_probabilities` (:class:`dict` of :class:`str`::class:`float`) -
          the probability distribution of the final number of main stem leaves.
          The keys are the final numbers of main stem leaves, and the values are
          the probabilities distribution.

        - `ears_density` (:class:`int`) - the number of ears per square meter, or None if no regression occured

        - `GL_number` (:class:`dict` of :class:`float`::class:`float`) - the GL decimal numbers measured at
          several thermal times (including the senescence end). The keys are the
          thermal times, and the values are the GL decimal numbers.

        - `delais_TT_stop_del_axis` (:class:`int`) - This variable represents the time in
          thermal time between an axis stop growing and its disappearance (it
          concerns only the axes that do not regress and which do not produce any
          cob).

        - `TT_hs_break` (:class:`float`) - the thermal time when the rate of Haun Stage
          is changing. *TT_hs_break* equal `None` (the default) means that the phyllochron is constant.

        - `inner_params` (:class:`dict`) - the values of the inner parameters used
          for the construction of the input tables. These parameters are the same
          as the ones defined in the module :mod:`params <openalea.adel.plantgen.params>`.
          *inner_params* is NOT mandatory: if not all inner parameters are documented
          in *inner_params*, then we use the default values defined in :mod:`params <openalea.adel.plantgen.params>`
          for the inner parameters which are missing.

        - axeT_user (:class:`pandas.DataFrame`): a table similar to the axeT_tmp that allows forcing which axis should be reconstructed.

        - TT_regression_start_user (:class: `float`) : thermal time at wich regression start on most frequent MS. If set to none, TT_regression_start is computed by pgen from start of MS elongation date of most frequent MS

        - TT_t1_user (:class: `float`) : thermal time at which n1 Green leaves are observed on most frequent MS. If set to none, TT_t1_user is computed by pgen from start of MS elongation date of most frequent MS

    :Returns:
        Return :ref:`axeT <axeT>`, :ref:`dimT <dimT>`,
        :ref:`phenT <phenT>`, :ref:`phenT_abs <phenT_abs>`,
        :ref:`dynT <dynT>`, :ref:`phenT_first <phenT_first>`, :ref:`HS_GL_SSI_T <HS_GL_SSI_T>`,
        :ref:`tilleringT <tilleringT>`, :ref:`cardinalityT <cardinalityT>`, and
        a dictionary which stores the configuration used for the construction.

    :Returns Type:
        tuple

    .. seealso:: :mod:`openalea.adel.plantgen.plantgen_core`
                 :mod:`openalea.adel.plantgen.params`
                 :mod:`openalea.adel.plantgen.tools`

    """
    # save the name and the value of each argument
    config = locals()

    # update values defined in openalea.adel.plantgen.params from values in inner_params
    attribute_names = set(dir(params))
    attribute_names.intersection_update(list(inner_params.keys()))
    params.__dict__.update(
        dict(
            [
                (key, value)
                for key, value in list(inner_params.items())
                if key in attribute_names
            ]
        )
    )

    if sum(MS_leaves_number_probabilities.values()) != 1.0:
        raise tools.InputError(
            "the sum of the probabilities defined in MS_leaves_number_probabilities is not equal to 1.0"
        )
    possible_axes = set(
        [
            id_axis
            for (id_axis, probability) in decide_child_axis_probabilities.items()
            if probability != 0.0
        ]
    )

    possible_MS_N_phytomer_potential = set(
        [
            int(MS_N_phytomer_potential)
            for (
                MS_N_phytomer_potential,
                probability,
            ) in MS_leaves_number_probabilities.items()
            if probability != 0.0
        ]
    )

    # check plants_number, decide_child_axis_probabilities, plants_density and ears_density validity
    decide_child_cohort_probabilities = (
        tools.calculate_decide_child_cohort_probabilities(
            decide_child_axis_probabilities
        )
    )
    (theoretical_cohort_cardinalities, theoretical_axis_cardinalities) = (
        tools.calculate_theoretical_cardinalities(
            plants_number,
            decide_child_cohort_probabilities,
            decide_child_axis_probabilities,
            params.FIRST_CHILD_DELAY,
            params.EMERGENCE_PROBABILITY_REDUCTION_FACTOR,
        )
    )
    theoretical_cardinalities_sum = sum(theoretical_cohort_cardinalities.values())
    if ears_density is not None:
        number_of_ears = plants_number * ears_density / float(plants_density)
    else:
        number_of_ears = None
    if number_of_ears is not None:
        if number_of_ears < plants_number:
            raise tools.InputError(
                "the number of ears (%s) is lesser than plants_number (%s). \
    The number of ears (%s) is calculated from plants_number (%s), ears_density (%s) \
    and plants_density (%s)."
                % (
                    number_of_ears,
                    plants_number,
                    number_of_ears,
                    plants_number,
                    ears_density,
                    plants_density,
                )
            )

        if number_of_ears > theoretical_cardinalities_sum:
            raise tools.InputError(
                "the number of ears (%s) is greater than the theoretical number of axes (%s). \
    The number of ears (%s) is calculated from plants_number (%s), ears_density (%s) and plants_density (%s). \
    The theoretical number of axes (%s) is calculated from plants_number (%s), decide_child_cohort_probabilities \
    (%s) and params.FIRST_CHILD_DELAY (%s)"
                % (
                    number_of_ears,
                    theoretical_cardinalities_sum,
                    number_of_ears,
                    plants_number,
                    ears_density,
                    plants_density,
                    theoretical_cardinalities_sum,
                    plants_number,
                    decide_child_cohort_probabilities,
                    params.FIRST_CHILD_DELAY,
                )
            )

    available_axes_warning_message = (
        "the probabilities defined in decide_child_axis_probabilities (%s) and \
the axes documented by the user (%s) in %s indicate that some of the possible axes (%s) \
are not documented by the user. After the generation of the axes, if not all generated axes are documented by the user, \
then this will lead to an error."
    )

    available_MS_N_phytomer_potential_warning_message = (
        "the probabilities defined in MS_leaves_number_probabilities (%s) and \
the N_phytomer_potential of the MS documented by the user (%s) in %s indicate that some of the N_phytomer_potential of the MS (%s) \
are not documented by the user. After the generation of the phytomers of the MS, if not all generated phytomers \
of the MS are documented by the user, then this will lead to an error."
    )

    # check the consistency of decide_child_axis_probabilities and params.MS_HS_AT_TILLER_EMERGENCE
    available_axes = list(params.MS_HS_AT_TILLER_EMERGENCE.keys())
    if not possible_axes.issubset(set(available_axes)):
        warnings.warn(
            available_axes_warning_message
            % (
                decide_child_axis_probabilities,
                available_axes,
                "params.MS_HS_AT_TILLER_EMERGENCE",
                list(possible_axes),
            ),
            tools.InputWarning,
        )

    # calculate dynT_user completeness and check its validity
    # update TT_col_* names with new names
    if "N_phytomer_potential" in dynT_user.columns:
        dynT_user_completeness = plantgen_core.DataCompleteness.FULL
        expected_dynT_user_columns = [
            "id_axis",
            "N_phytomer_potential",
            "a_cohort",
            "TT_hs_0",
            "TT_flag_ligulation",
            "n0",
            "n1",
            "n2",
        ]
        if set(dynT_user.columns.tolist()) != set(expected_dynT_user_columns):
            raise tools.InputError(
                "dynT_user does not have the columns: %s"
                % ", ".join(expected_dynT_user_columns)
            )
        grouped = dynT_user.groupby(["id_axis", "N_phytomer_potential"])
        if len(grouped.groups) != dynT_user.index.size:
            raise tools.InputError(
                "dynT_user contains duplicated (id_axis, N_phytomer_potential) pair(s)"
            )
        available_axes = dynT_user["id_axis"].tolist()
        if not possible_axes.issubset(set(available_axes)):
            warnings.warn(
                available_axes_warning_message
                % (
                    decide_child_axis_probabilities,
                    available_axes,
                    "dynT_user",
                    list(possible_axes),
                ),
                tools.InputWarning,
            )
        available_MS_N_phytomer_potential = set(
            dynT_user[dynT_user["id_axis"] == "MS"]["N_phytomer_potential"].tolist()
        )
        if not possible_MS_N_phytomer_potential.issubset(
            available_MS_N_phytomer_potential
        ):
            warnings.warn(
                available_MS_N_phytomer_potential_warning_message
                % (
                    MS_leaves_number_probabilities,
                    list(available_MS_N_phytomer_potential),
                    "dynT_user",
                    list(possible_MS_N_phytomer_potential),
                ),
                tools.InputWarning,
            )

    elif dynT_user.count().max() == dynT_user.count().min() == dynT_user.index.size:
        dynT_user_completeness = plantgen_core.DataCompleteness.SHORT
        expected_dynT_user_columns = [
            "id_axis",
            "a_cohort",
            "TT_hs_0",
            "TT_flag_ligulation",
            "n0",
            "n1",
            "n2",
        ]
        if set(dynT_user.columns.tolist()) != set(expected_dynT_user_columns):
            raise tools.InputError(
                "dynT_user does not have the columns: %s"
                % ", ".join(expected_dynT_user_columns)
            )
        if dynT_user["id_axis"].unique().size != dynT_user["id_axis"].size:
            raise tools.InputError("dynT_user contains duplicated id_axis")
        available_axes = dynT_user["id_axis"].tolist()
        if not possible_axes.issubset(set(available_axes)):
            warnings.warn(
                available_axes_warning_message
                % (
                    decide_child_axis_probabilities,
                    available_axes,
                    "dynT_user",
                    list(possible_axes),
                ),
                tools.InputWarning,
            )

    else:
        dynT_user_completeness = plantgen_core.DataCompleteness.MIN
        expected_dynT_user_columns = [
            "id_axis",
            "a_cohort",
            "TT_hs_0",
            "TT_flag_ligulation",
            "n0",
            "n1",
            "n2",
        ]
        if set(dynT_user.columns.tolist()) != set(expected_dynT_user_columns):
            raise tools.InputError(
                "dynT_user does not have the columns: %s"
                % ", ".join(expected_dynT_user_columns)
            )
        available_axes = dynT_user["id_axis"].tolist()
        if not possible_axes.issubset(set(available_axes)):
            warnings.warn(
                available_axes_warning_message
                % (
                    decide_child_axis_probabilities,
                    available_axes,
                    "dynT_user",
                    list(possible_axes),
                ),
                tools.InputWarning,
            )

    # calculate dimT_user completeness and check its validity
    if "N_phytomer_potential" in dimT_user.columns:
        dimT_user_completeness = plantgen_core.DataCompleteness.FULL
        expected_dimT_user_columns = [
            "id_axis",
            "N_phytomer_potential",
            "index_phytomer",
            "L_blade",
            "W_blade",
            "L_sheath",
            "W_sheath",
            "L_internode",
            "W_internode",
        ]
        if set(dimT_user.columns.tolist()) != set(expected_dimT_user_columns):
            raise tools.InputError(
                "dimT_user does not have the columns: %s"
                % ", ".join(expected_dimT_user_columns)
            )
        grouped = dimT_user.groupby(
            ["id_axis", "N_phytomer_potential", "index_phytomer"]
        )
        if len(grouped.groups) != dimT_user.index.size:
            raise tools.InputError(
                "dimT_user contains duplicated (id_axis, N_phytomer_potential, index_phytomer) triplet(s)"
            )
        available_axes = dimT_user["id_axis"].tolist()
        if not possible_axes.issubset(set(available_axes)):
            warnings.warn(
                available_axes_warning_message
                % (
                    decide_child_axis_probabilities,
                    available_axes,
                    "dimT_user",
                    list(possible_axes),
                ),
                tools.InputWarning,
            )
        available_MS_N_phytomer_potential = set(
            dimT_user[dimT_user["id_axis"] == "MS"]["N_phytomer_potential"].tolist()
        )
        if not possible_MS_N_phytomer_potential.issubset(
            available_MS_N_phytomer_potential
        ):
            warnings.warn(
                available_MS_N_phytomer_potential_warning_message
                % (
                    MS_leaves_number_probabilities,
                    list(available_MS_N_phytomer_potential),
                    "dimT_user",
                    list(possible_MS_N_phytomer_potential),
                ),
                tools.InputWarning,
            )

    elif "id_axis" in dimT_user.columns:
        dimT_user_completeness = plantgen_core.DataCompleteness.SHORT
        expected_dimT_user_columns = [
            "id_axis",
            "index_phytomer",
            "L_blade",
            "W_blade",
            "L_sheath",
            "W_sheath",
            "L_internode",
            "W_internode",
        ]
        if set(dimT_user.columns.tolist()) != set(expected_dimT_user_columns):
            raise tools.InputError(
                "dimT_user does not have the columns: %s"
                % ", ".join(expected_dimT_user_columns)
            )
        grouped = dimT_user.groupby(["id_axis", "index_phytomer"])
        if len(grouped.groups) != dimT_user.index.size:
            raise tools.InputError(
                "dimT_user contains duplicated (id_axis, index_phytomer) pair(s)"
            )
        available_axes = dimT_user["id_axis"].tolist()
        if not possible_axes.issubset(set(available_axes)):
            warnings.warn(
                available_axes_warning_message
                % (
                    decide_child_axis_probabilities,
                    available_axes,
                    "dimT_user",
                    list(possible_axes),
                ),
                tools.InputWarning,
            )
        max_available_MS_N_phytomer_potential = dimT_user[dimT_user["id_axis"] == "MS"][
            "index_phytomer"
        ].max()
        if (
            max(possible_MS_N_phytomer_potential)
            > max_available_MS_N_phytomer_potential
        ):
            warnings.warn(
                available_MS_N_phytomer_potential_warning_message
                % (
                    MS_leaves_number_probabilities,
                    ", ".join([str(max_available_MS_N_phytomer_potential)]),
                    "dimT_user",
                    ", ".join([str(max(possible_MS_N_phytomer_potential))]),
                ),
                tools.InputWarning,
            )

    else:
        dimT_user_completeness = plantgen_core.DataCompleteness.MIN
        expected_dimT_user_columns = [
            "index_phytomer",
            "L_blade",
            "W_blade",
            "L_sheath",
            "W_sheath",
            "L_internode",
            "W_internode",
        ]
        if set(dimT_user.columns.tolist()) != set(expected_dimT_user_columns):
            raise tools.InputError(
                "dimT_user does not have the columns: %s"
                % ", ".join(expected_dimT_user_columns)
            )
        if (
            dimT_user["index_phytomer"].unique().size
            != dimT_user["index_phytomer"].size
        ):
            raise tools.InputError("dimT_user contains duplicated index_phytomer")
        max_available_MS_N_phytomer_potential = dimT_user["index_phytomer"].max()
        if (
            max(possible_MS_N_phytomer_potential)
            > max_available_MS_N_phytomer_potential
        ):
            warnings.warn(
                available_MS_N_phytomer_potential_warning_message
                % (
                    MS_leaves_number_probabilities,
                    ", ".join([str(max_available_MS_N_phytomer_potential)]),
                    "dimT_user",
                    ", ".join([str(max(possible_MS_N_phytomer_potential))]),
                ),
                tools.InputWarning,
            )

    if TT_hs_break is None:
        TT_hs_break = np.nan

    # initialize the axes
    cardinalityT = plantgen_core.init_axes(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        theoretical_cohort_cardinalities,
        theoretical_axis_cardinalities,
        axeT_user=axeT_user,
    )

    # define phenology functions
    dynT_, decimal_elongated_internode_number = plantgen_core.phenology_functions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        dynT_user_completeness,
        dimT_user_completeness,
        TT_hs_break,
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )

    # construct plants structure
    axeT_, tilleringT, phenT_first = plantgen_core.plants_structure(
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
        axeT_user=axeT_user,
        TT_regression_start_user=TT_regression_start_user,
        TT_t1_user=TT_t1_user,
    )

    # compute organs dimensions
    dimT_ = plantgen_core.organs_dimensions(
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
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )
    # calculate the phenology of the axes
    phenT_, phenT_abs, HS_GL_SSI_T = plantgen_core.axes_phenology(
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
        axeT_user=axeT_user,
        TT_t1_user=TT_t1_user,
    )

    return (
        axeT_,
        dimT_,
        phenT_,
        phenT_abs,
        dynT_,
        phenT_first,
        HS_GL_SSI_T,
        tilleringT,
        cardinalityT,
        config,
    )


def read_plantgen_inputs(inputs_filepath, dynT_user_filepath, dimT_user_filepath):
    """
    Import the Python module at *inputs_filepath*, and return the args expected by
    :func:`gen_adel_input_data`.

    :Parameters:

        - `inputs_filepath` (:class:`str`) - the file path of the Python module
          which contains the inputs of :func:`gen_adel_input_data`.

        - `dynT_user_filepath` (:class:`str`) - the file path of
          the leaf dynamic parameters set by the user.

        - `dimT_user_filepath` (:class:`str`) - the file path of
          the dimensions of the organs set by the user.

    :Returns:
        Return the inputs of :func:`gen_adel_input_data`.

    :Returns Type:
        tuple

    """
    import imp

    inputs = imp.load_source("inputs", inputs_filepath)

    dynT_user = pd.read_csv(dynT_user_filepath)
    dimT_user = pd.read_csv(dimT_user_filepath)
    plants_number = inputs.plants_number
    plants_density = inputs.plants_density
    decide_child_axis_probabilities = inputs.decide_child_axis_probabilities
    MS_leaves_number_probabilities = inputs.MS_leaves_number_probabilities
    ears_density = inputs.ears_density
    GL_number = inputs.GL_number
    delais_TT_stop_del_axis = inputs.delais_TT_stop_del_axis

    try:
        TT_hs_break = inputs.TT_hs_break
    except:
        TT_hs_break = None

    try:
        inner_params = inputs.inner_params
    except:
        inner_params = {}

    return (
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
    )


def plantgen2adel(axeT_, dimT_, phenT_):
    """
    Format the dataframes generated by :mod:`gen_adel_input_data <openalea.adel.plantgen.plantgen_interface.gen_adel_input_data>`
    to dataframes compatible with adel. This function is a temporary patch, waiting
    for Adel to be updated.

    :Parameters:

        - `axeT_` (:class:`pandas.DataFrame`) - the :ref:`axeT <axeT>` dataframe.
        - `dimT_` (:class:`pandas.DataFrame`) - the :ref:`dimT <dimT>` dataframe.
        - `phenT_` (:class:`pandas.DataFrame`) - the :ref:`phenT <phenT>` dataframe.

    :Returns:
        Return :ref:`axeT <axeT>`, :ref:`dimT <dimT>` and :ref:`phenT <phenT>`
        in adel-like format.

    :Returns Type:
        tuple

    """

    # In dimT, duplicate the line for each axis which has only one leaf.
    id_dims_with_one_leaf = axeT_[axeT_["N_phytomer"] == 1]["id_dim"].unique()
    for id_dim in id_dims_with_one_leaf:
        idx = dimT_[dimT_["id_dim"] == id_dim].first_valid_index()
        new_line = dimT_.ix[idx:idx].copy()
        new_line["index_rel_phytomer"] = 0.5
        dimT_ = pd.concat([dimT_[:idx], new_line, dimT_[idx:]], ignore_index=True)

    return axeT_, dimT_, phenT_
