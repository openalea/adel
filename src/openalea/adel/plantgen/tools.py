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
Generic tools used in the :mod:`openalea.adel.plantgen` package. These routines can
also be used by other packages.

Authors: M. Abichou, B. Andrieu, C. Chambon
"""

import random

import math

import numpy as np
import pandas as pd
from scipy.optimize import leastsq


def decide_child_cohorts(
    decide_child_cohort_probabilities,
    first_child_delay,
    emergence_probability_reduction_factor,
    parent_cohort_index=None,
    parent_cohort_position=None,
):
    """
    Decide (recursively) of the child cohorts actually produced by a parent cohort,
    according to the *decide_child_cohort_probabilities* and the *parent_cohort_index*.
    The main stem always exists.

    :Parameters:

        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities of the
          child cohorts.
        - `first_child_delay` (:class:`int`) - the delay between the parent cohort and
          the first child cohort. This delay is expressed in number of cohorts.
        - `emergence_probability_reduction_factor` (:class:`float`) - The reduction factor
          of the emergence probability of secondary tiller compared to primary one.
        - `parent_cohort_index` (:class:`int`) - the index of the parent cohort.
          ``None`` (the default) means that there isn't any parent cohort.
        - `parent_cohort_position` (:class:`str`) - the position of the parent cohort.
          ``None`` (the default) means that there isn't any parent cohort.

    :Returns:
        The indices of the child cohorts and their positions in the tree.

    :Returns Type:
        list of tuples

    """
    child_cohorts = []
    if parent_cohort_index is None:
        first_possible_cohort_number = 1
    else:
        first_possible_cohort_number = parent_cohort_index + first_child_delay
    if first_possible_cohort_number == 1:
        # The main stem always exists: add it.
        cohort_position = "MS"
        child_cohorts.append((first_possible_cohort_number, "MS"))
        child_cohorts.extend(
            decide_child_cohorts(
                decide_child_cohort_probabilities,
                first_child_delay,
                emergence_probability_reduction_factor,
                first_possible_cohort_number,
                cohort_position,
            )
        )
    else:
        # Find the children of the secondary stem.
        emergence_probability_reduction_coefficient = (
            1.0 - emergence_probability_reduction_factor
        )
        for cohort_id, cohort_probability in decide_child_cohort_probabilities.items():
            if cohort_id >= first_possible_cohort_number:
                reducted_cohort_probability = cohort_probability
                if parent_cohort_position != "MS":
                    reducted_cohort_probability *= (
                        emergence_probability_reduction_coefficient  # apply a reduction
                    )
                if reducted_cohort_probability >= random.random():
                    child_cohort_position = (
                        cohort_id - parent_cohort_index - first_child_delay
                    )
                    if parent_cohort_position == "MS":
                        cohort_position = "T%s" % child_cohort_position
                    else:
                        cohort_position = "%s.%s" % (
                            parent_cohort_position,
                            child_cohort_position,
                        )
                    child_cohorts.append((cohort_id, cohort_position))
                    child_cohorts.extend(
                        decide_child_cohorts(
                            decide_child_cohort_probabilities,
                            first_child_delay,
                            emergence_probability_reduction_factor,
                            cohort_id,
                            cohort_position,
                        )
                    )
    return child_cohorts


def calculate_MS_final_leaves_number(MS_leaves_number_probabilities):
    """
    Calculate the final number of leaves of a main stem. This is done by randomly
    drawing a number in a probability distribution. Uses the probabilities
    of the main stem leaves number.

    :Parameters:

        - `MS_leaves_number_probabilities` (:class:`dict`) - the probabilities
          of the main stem leaves number.

    :Returns:
        The final number of leaves of the main stem.

    :Returns Type:
        :class:`float`

    """
    random_value = random.random()
    probabilities_sum = 0.0
    MS_final_leaves_number = None
    for leaves_number_str, leaves_probability in MS_leaves_number_probabilities.items():
        probabilities_sum += leaves_probability
        if random_value <= probabilities_sum:
            MS_final_leaves_number = float(leaves_number_str)
            break
    return MS_final_leaves_number


def calculate_tiller_final_leaves_number(
    MS_final_leaves_number, cohort_number, secondary_stem_leaves_number_coefficients
):
    """
    Calculate the final number of leaves of a tiller.
    Uses the final number of leaves of  the main stem, the index of the cohort to
    which belongs the tiller, and specific coefficients.

    :Parameters:

        - `MS_final_leaves_number` (:class:`float`) - the final number of leaves of the
          main stem.
        - `cohort_number` (:class:`int`) - the index of cohort.
        - `secondary_stem_leaves_number_coefficients` (:class:`dict`) - The coefficients
          a_1 and a_2 to calculate the final number of leaves on tillers from
          the final number of leaves on main stem. Calculation is done as follow::

              tiller_final_leaves_number
                  = a_1 * MS_final_leaves_number - a_2 * cohort_number

    :Returns:
        The final number of leaves of a tiller.

    :Returns Type:
        :class:`float`

    """
    a_1 = secondary_stem_leaves_number_coefficients["a_1"]
    a_2 = secondary_stem_leaves_number_coefficients["a_2"]
    return a_1 * MS_final_leaves_number - a_2 * cohort_number


def decide_time_of_death(
    max_axes_number,
    number_of_ears,
    TT_regression_start,
    TT_regression_end,
    TT_em_phytomer1_df,
):
    """
    Decide the thermal times (relative to canopy emergence) when the axes stop
    growing. Uses an exponential function which describes the decay of the global population.

    :Parameters:

        - `max_axes_number` (:class:`int`) - the maximum number of existing axes.

        - `number_of_ears` (:class:`int`) - the number of ears.

        - `TT_regression_start` (:class:`float`) - thermal time at which the regression starts, i.e. when the Haun Stage
            of the most frequent main stem is equal to (N_phytomer_potential - 5).

        - `TT_regression_end` (:class:`float`) - the thermal time at which the regression ends, i.e. when the haun stage
            of the most frequent main stem is equal to flag leaf number.

        - `TT_em_phytomer1_df` (:class:`pandas.DataFrame`) - A dataframe which contains, for each plant and each tiller, the thermal time
            (relative to canopy appearance) of tip appearance of the first true leaf (not coleoptile nor prophyll).

    :Returns:
        the thermal times (relative to canopy emergence) when the axes stop growing.

    :Returns Type:
        :class:`list`

    .. warning::

        * *number_of_ears*, *max_axes_number*, *TT_regression_start* and *TT_regression_end*
          must be positive or null.
        * *TT_regression_start* must be smaller (or equal) than *TT_regression_end*.
        * *number_of_ears* must be smaller (or equal) than *max_axes_number*.

    """
    if number_of_ears is None:  # no regression
        return [np.nan] * len(TT_em_phytomer1_df)
    else:
        if max_axes_number < 0:
            raise InputError("max_axes_number negative")
        if number_of_ears < 0:
            raise InputError("number_of_ears negative")
        if TT_regression_start < 0:
            raise InputError("TT_regression_start negative")
        if TT_regression_end < 0:
            raise InputError("TT_regression_end negative")

        if TT_regression_start > TT_regression_end:
            raise InputError("TT_regression_start greater than TT_regression_end")

        if number_of_ears > max_axes_number:
            raise InputError("number_of_ears greater than max_axes_number")

        def calculate_number_of_active_axes(
            tt, max_axes_number, number_of_ears, TT_regression_start, TT_regression_end
        ):
            return (max_axes_number - number_of_ears) * math.exp(
                -2.59861720216721
                * (tt - TT_regression_start)
                / (
                    1.65412664908155 * (TT_regression_end - TT_regression_start)
                    - (tt - TT_regression_start)
                )
            ) + number_of_ears  ## les constantes sont a mettre dans params

        TT_stop_axis_df = TT_em_phytomer1_df.assign(
            TT_stop_axis=pd.Series(index=TT_em_phytomer1_df.index)
        )
        number_of_remaining_axes = max_axes_number

        for tt in range(int(TT_regression_start), int(TT_regression_end) + 1):
            number_of_active_axes = int(
                calculate_number_of_active_axes(
                    tt,
                    max_axes_number,
                    number_of_ears,
                    TT_regression_start,
                    TT_regression_end,
                )
            )
            number_of_axes_to_delete = number_of_remaining_axes - number_of_active_axes
            # we consider only the lines where TT_stop_axis is still NA
            TT_stop_axis_NA_df = TT_stop_axis_df.loc[
                TT_stop_axis_df.TT_stop_axis.isnull()
            ]
            plants_to_use_for_random_draw = TT_stop_axis_NA_df.id_plt.unique().tolist()
            while number_of_axes_to_delete > 0:
                id_plt = random.choice(plants_to_use_for_random_draw)
                TT_stop_axis_NA_df_plant_group = TT_stop_axis_NA_df.loc[
                    TT_stop_axis_NA_df.id_plt == id_plt
                ]
                index_of_yougest_tiller = (
                    TT_stop_axis_NA_df_plant_group.TT_em_phytomer1.idxmax()
                )
                TT_stop_axis_df.loc[index_of_yougest_tiller, "TT_stop_axis"] = tt
                plants_to_use_for_random_draw.remove(id_plt)
                number_of_axes_to_delete -= 1
                number_of_remaining_axes -= 1
            if number_of_remaining_axes == 0:
                break
        return TT_stop_axis_df.TT_stop_axis.tolist()


def fit_poly(x_meas_array, y_meas_array, fixed_coefs, a_starting_estimate):
    """
    Calculate the best-fit parameter *a*, where *a* is the coefficient of highest
    degree of a polynomial. The other polynomial coefficients are supposed to be
    known and are given in *fixed_coefs*.
    We first define a function to compute the residuals. Then we use the least-squares
    fit routine of scipy to find the best-fit parameter *a*, selecting *a_starting_estimate*
    as starting position and using (*x_meas_array*, *y_meas_array*) as the measured
    data to fit. Finally, we calculate the *RMSE* to check the validity of the fit.

    .. seealso:: :func:`scipy.optimize.leastsq`

    :Parameters:

        - `x_meas_array` (:class:`np.ndarray`) - the x-coordinates. These data are
          measured.
        - `y_meas_array` (:class:`np.ndarray`) - the y-coordinates. These data are
          measured.
        - `fixed_coefs` (:class:`list`) - the other coefficients of the polynomial to fit
          (*x_meas_array*, *y_meas_array*) to.
          These coefficients are not fitted. They are given from highest degree
          to lowest degree ("descending powers").
        - `a_starting_estimate` (float) - the starting estimate for the minimization.

    :Returns:
        the best-fit coefficient *a* and the *RMSE* of the fit.

    :Returns Type:
        :class:`tuple` of :class:`float`

    """

    def residuals(p, y, x):
        (a,) = p
        err = y - peval(x, a)
        return err

    def peval(x, a):
        return np.poly1d([a] + fixed_coefs)(x)

    p, cov, infodict, mesg, ier = leastsq(
        residuals,
        [a_starting_estimate],
        args=(y_meas_array, x_meas_array),
        full_output=1,
    )
    # RMSE_gl
    chisq = (infodict["fvec"] ** 2).sum()
    dof = len(x_meas_array) - 1  # dof is degrees of freedom
    rmse = np.sqrt(chisq / dof)
    return p[0], rmse


def calculate_theoretical_cardinalities(
    plants_number,
    decide_child_cohort_probabilities,
    decide_child_axis_probabilities,
    first_child_delay,
    emergence_probability_reduction_factor,
):
    """
    Calculate the theoretical cardinality of each simulated cohort and each
    simulated axis.

    :Parameters:

        - `plants_number` (:class:`int`) - the number of plants.
        - `decide_child_cohort_probabilities` (:class:`dict`) - the probabilities
          of the child cohorts.
        - `decide_child_axis_probabilities` (:class:`dict`) - the probabilities
          of the child axes.
        - `first_child_delay` (:class:`int`) - The delay between
          a parent axis and its first possible child axis. This delay is
          expressed in number of cohorts.
        - `emergence_probability_reduction_factor` (:class:`float`) - The reduction factor
          of the emergence probability of secondary tiller compared to primary one.

    :Returns:
        a 2-tuple of dictionaries: the first dictionary contains the theoretical
        cardinality of each cohort, the second dictionary contains the theoretical
        cardinality of each axis.

    :Returns Type:
        :class:`tuple`

    """
    child_cohort_probabilities_ceiled = dict(
        list(
            zip(
                list(decide_child_cohort_probabilities.keys()),
                np.ceil(list(decide_child_cohort_probabilities.values())),
            )
        )
    )
    all_child_cohorts = set()
    for i in range(plants_number):
        child_cohorts = decide_child_cohorts(
            child_cohort_probabilities_ceiled, first_child_delay, 0.0
        )
        all_child_cohorts.update(child_cohorts)

    axis_to_cohort_mapping = dict(
        [(cohort_axis[1], cohort_axis[0]) for cohort_axis in all_child_cohorts]
    )

    id_cohort_list = sorted([1] + list(decide_child_cohort_probabilities.keys()))
    theoretical_cohort_cardinalities = dict.fromkeys(id_cohort_list, 0.0)
    id_axis_list = sorted(["MS"] + list(decide_child_axis_probabilities.keys()))
    id_cohort_id_axis_tuples = list(zip(id_cohort_list, id_axis_list))
    theoretical_axis_cardinalities = dict.fromkeys(id_cohort_id_axis_tuples, 0.0)
    emergence_probability_reduction_coefficient = (
        1.0 - emergence_probability_reduction_factor
    )
    for id_cohort, id_axis in all_child_cohorts:
        if id_cohort == 1:  # main stem
            axis_probability = 1.0
        else:  # tillers
            axis_probability = decide_child_cohort_probabilities[id_cohort]
            current_id_axis = id_axis
            while "." in current_id_axis:
                current_parent_axis = current_id_axis.rsplit(".", 1)[0]
                axis_probability *= (
                    decide_child_cohort_probabilities[
                        axis_to_cohort_mapping[current_parent_axis]
                    ]
                    * emergence_probability_reduction_coefficient
                )
                current_id_axis = current_parent_axis

        number_of_axes = axis_probability * plants_number
        theoretical_cohort_cardinalities[id_cohort] += number_of_axes
        theoretical_axis_cardinalities[(id_cohort, id_axis)] = number_of_axes

    return (theoretical_cohort_cardinalities, theoretical_axis_cardinalities)


def calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities):
    """
    For each primary tiller in *decide_child_axis_probabilities*, calculate the corresponding
    cohort number, and return a dictionary which keys are cohort number and values remain
    the same.

    :Parameters:

        - `decide_child_axis_probabilities` (:class:`dict`) - the probability for each
          primary tiller to have a child. Keys are the botanical positions (e.g. "T1", "T2",...),
          values are the probabilities (float).

    :Returns:
        the probability for each cohort to have a child. Keys are the indexes of the cohorts
        (e.g. 3, 4,...), values are the probabilities (float).

    :Returns Type:
        :class:`dict`

    :Examples:

        >>> decide_child_axis_probabilities = {'T0': 0.0, 'T1': 0.900, 'T2': 0.983, 'T3': 0.817, 'T4': 0.117}
        >>> calculate_decide_child_cohort_probabilities(decide_child_axis_probabilities)
        {3: 0.0, 4: 0.900, 5: 0.983, 6: 0.817, 7: 0.117}

    """
    id_axis_array = np.array(list(decide_child_axis_probabilities.keys()))
    id_cohort_array = np.char.lstrip(id_axis_array, "T").astype(int) + 3
    decide_child_cohort_probabilities = dict(
        list(zip(id_cohort_array, list(decide_child_axis_probabilities.values())))
    )
    return decide_child_cohort_probabilities


def get_primary_axis(id_axis, first_child_delay):
    """
    Calculate the primary axis of *id_axis*.

    :Parameters:

        - `id_axis` (:class:`str`) - the botanical position of the axis.
        - `first_child_delay` (:class:`int`) - the delay between the axis and its
          first child.

    :Returns:
        the primary axis of *id_axis*.

    :Returns Type:
        :class:`str`

    :Examples:

        >>> get_primary_axis('T1.0', 2)
        'T3'
        >>> get_primary_axis('T1.0.0', 2)
        'T5'
        >>> get_primary_axis('T5', 2)
        'T5'
        >>> get_primary_axis('MS', 2)
        'MS'

    """
    if id_axis == "MS":
        primary_id_axis = id_axis
    else:
        id_axis = id_axis[1:]
        while "." in id_axis:
            id_axis_split = id_axis.rsplit(".", 2)
            last_pos = int(id_axis_split.pop())
            last_but_one_pos = int(id_axis_split.pop())
            new_last_pos = last_but_one_pos + last_pos + first_child_delay
            id_axis = ".".join(id_axis_split + [str(new_last_pos)])
        primary_id_axis = "T" + id_axis
    return primary_id_axis


def get_real_roots(poly):
    """
    Get the real roots of polynomial *poly*.

    :Parameters:

        - `poly` (:class:`numpy.lib.polynomial.poly1d`) - a one-dimensional polynomial.

    :Returns:
        the real roots of *poly*

    :Returns Type:
        :class:`numpy.array`

    :Examples:

        >>> p = numpy.poly1d([1, 2, 1])
        array([-1., -1.])
        >>> p = numpy.poly1d([1, 2, 3])
        array([], dtype=float64)

    """
    roots_array = poly.r
    real_roots = [x for x in roots_array if x.imag == 0.0]
    real_roots = [x.real for x in real_roots]
    return np.array(real_roots)


def find_lines_intersection(line1, line2):
    """
    Find the intersection of two lines.
    Raise an exception if no intersection.

    :Parameters:

        - `line1` (:class:`tuple` of 2 :class:`tuple` of :class:`float`) - the coordinates of the two points which define the first line.
        - `line2` (:class:`tuple` of 2 :class:`tuple` of :class:`float`) - the coordinates of the two points which define the second line.

    :Returns:
        the coordinates of the intersection point.

    :Returns Type:
        :class:`tuple` of :class:`float`

    :Examples:

        >>> find_lines_intersection(((0.5, 0.5), (1.5, 0.5)), ((0.5, 0.5), (1.5, 0.5)))
        (1.0, 0.5)

    .. codeauthor:: from http://stackoverflow.com/a/20679579

    """

    def compute_line_equation_coefficients(point1, point2):
        # Compute coefficients a, b and c of line equation by two points provided.
        a = point1[1] - point2[1]
        b = point2[0] - point1[0]
        c = -(point1[0] * point2[1] - point2[0] * point1[1])
        return a, b, c

    a1, b1, c1 = compute_line_equation_coefficients(line1[0], line1[1])
    a2, b2, c2 = compute_line_equation_coefficients(line2[0], line2[1])

    main_determinant = a1 * b2 - b1 * a2
    x_determinant = c1 * b2 - b1 * c2
    y_determinant = a1 * c2 - c1 * a2

    if main_determinant == 0:
        raise Exception("Lines do not intersect.")

    intersection_point = (
        x_determinant / main_determinant,
        y_determinant / main_determinant,
    )
    return intersection_point


class InputError(Exception):
    """Exception raised when an invalid input is detected."""

    pass


class InputWarning(UserWarning):
    """Warning issued when an input is dubious and may lead to an error."""

    pass
