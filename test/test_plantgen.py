import random

from openalea.adel.plantgen import params, tools, plantgen_interface, plantgen_core
import numpy as np
import pandas
from pathlib import Path as path
import tempfile


random.seed(1234)

initial_random_state = random.getstate()

plants_number = 5
plants_density = 12
decide_child_axis_probabilities = {
    "T0": 0.0,
    "T1": 0.900,
    "T2": 0.967,
    "T3": 0.817,
    "T4": 0.083,
}
decide_child_cohort_probabilities = tools.calculate_decide_child_cohort_probabilities(
    decide_child_axis_probabilities
)
MS_leaves_number_probabilities = {
    "10": 0.145,
    "11": 0.818,
    "12": 0.037,
    "13": 0.0,
    "14": 0.0,
}
ears_density = 25
GL_number = {1117.0: 5.6, 1212.1: 5.4, 1368.7: 4.9, 1686.8: 2.4, 1880.0: 0.0}
delais_TT_stop_del_axis = 600
TT_hs_break = np.NaN  # linear mode
TT_flag_ligulation = {
    "MS": 1078.0,
    "T1": 1148.0,
    "T2": 1158.0,
    "T3": 1168.0,
    "T4": 1178.0,
}
number_of_ears = plants_number * ears_density / float(plants_density)

expected_results_dir = path(__file__).parent.resolve() / path("data/test_plantgen")
default_expected_results_dir = expected_results_dir.joinpath("default")
min_min_expected_results_dir = expected_results_dir.joinpath("min_min")
short_short_expected_results_dir = expected_results_dir.joinpath("short_short")
full_full_expected_results_dir = expected_results_dir.joinpath("full_full")

tmp_results_directory = path(tempfile.mkdtemp(suffix="_plantgen_results"))
default_results = tmp_results_directory.joinpath("default")
if not default_results.exists():
    default_results.mkdir()
relative_tolerance = 10e-3
absolute_tolerance = 10e-3

OUTPUTS_PRECISION = 10

FLOAT_FORMAT = "%.{}f".format(OUTPUTS_PRECISION)


def reinit_random_state():
    global initial_random_state
    random.setstate(initial_random_state)


def test_init_axes():
    reinit_random_state()
    (theoretical_cohort_cardinalities, theoretical_axis_cardinalities) = (
        tools.calculate_theoretical_cardinalities(
            plants_number,
            decide_child_cohort_probabilities,
            decide_child_axis_probabilities,
            params.FIRST_CHILD_DELAY,
            params.EMERGENCE_PROBABILITY_REDUCTION_FACTOR,
        )
    )
    cardinalityT = plantgen_core.init_axes(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        theoretical_cohort_cardinalities,
        theoretical_axis_cardinalities,
    )
    expected_cardinalityT = pandas.read_csv(
        default_expected_results_dir / "cardinalityT.csv"
    )
    test_table_filepath = default_results.joinpath("cardinalityT.csv")
    cardinalityT.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    assert (cardinalityT["id_axis"] == expected_cardinalityT["id_axis"]).all()
    cardinalityT = cardinalityT.drop("id_axis", axis=1)
    expected_cardinalityT = expected_cardinalityT.drop("id_axis", axis=1)
    np.testing.assert_allclose(
        cardinalityT.astype(float).values,
        expected_cardinalityT.astype(float).values,
        relative_tolerance,
        absolute_tolerance,
    )


def test_phenology_functions():
    reinit_random_state()
    dynT_user = pandas.read_csv(short_short_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(short_short_expected_results_dir / "dimT_user.csv")
    dynT_, decimal_elongated_internode_number = plantgen_core.phenology_functions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        plantgen_core.DataCompleteness.SHORT,
        plantgen_core.DataCompleteness.SHORT,
        TT_hs_break,
    )
    expected_dynT = pandas.read_csv(default_expected_results_dir / "dynT.csv")
    test_table_filepath = default_results.joinpath("dynT.csv")
    dynT_.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    assert len(dynT_) == len(expected_dynT)
    assert dynT_["id_axis"] == expected_dynT["id_axis"]
    dynT_ = dynT_.drop("id_axis", axis=1)
    expected_dynT = expected_dynT.drop(["id_axis"], axis=1)
    np.testing.assert_allclose(
        dynT_.values, expected_dynT.values, relative_tolerance, absolute_tolerance
    )


def test_plants_structure():
    reinit_random_state()
    dynT_user = pandas.read_csv(short_short_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(short_short_expected_results_dir / "dimT_user.csv")
    axeT_, tilleringT, phenT_first = plantgen_core.plants_structure(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        plantgen_core.DataCompleteness.SHORT,
        plantgen_core.DataCompleteness.SHORT,
        TT_hs_break,
        delais_TT_stop_del_axis,
        number_of_ears,
        plants_density,
        ears_density,
    )
    expected_axeT = pandas.read_csv(default_expected_results_dir / "axeT.csv")
    test_table_filepath = default_results.joinpath("axeT.csv")
    axeT_.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    assert (axeT_["id_axis"] == expected_axeT["id_axis"]).all()
    axeT_ = axeT_.drop(["id_axis"], axis=1)
    expected_axeT = expected_axeT.drop(["id_axis"], axis=1)
    np.testing.assert_allclose(
        axeT_.values, expected_axeT.values, relative_tolerance, absolute_tolerance
    )

    expected_tilleringT = pandas.read_csv(
        default_expected_results_dir / "tilleringT.csv"
    )
    test_table_filepath = default_results.joinpath("tilleringT.csv")
    tilleringT.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        tilleringT.values,
        expected_tilleringT.values,
        relative_tolerance,
        absolute_tolerance,
    )

    expected_phenT_first = pandas.read_csv(
        default_expected_results_dir / "phenT_first.csv"
    )
    test_table_filepath = default_results.joinpath("phenT_first.csv")
    phenT_first.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        phenT_first.values,
        expected_phenT_first.values,
        relative_tolerance,
        absolute_tolerance,
    )


def test_organs_dimensions():
    reinit_random_state()
    dynT_user = pandas.read_csv(short_short_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(short_short_expected_results_dir / "dimT_user.csv")

    dimT_ = plantgen_core.organs_dimensions(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        plantgen_core.DataCompleteness.SHORT,
        plantgen_core.DataCompleteness.SHORT,
        TT_hs_break,
        delais_TT_stop_del_axis,
        number_of_ears,
    )

    expected_dimT = pandas.read_csv(default_expected_results_dir / "dimT.csv")
    test_table_filepath = default_results.joinpath("dimT.csv")
    dimT_.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        dimT_.values.astype(float),
        expected_dimT.values,
        relative_tolerance,
        absolute_tolerance,
    )


def test_axes_phenology():
    reinit_random_state()
    dynT_user = pandas.read_csv(short_short_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(short_short_expected_results_dir / "dimT_user.csv")
    phenT_, phenT_abs, HS_GL_SSI_T = plantgen_core.axes_phenology(
        plants_number,
        decide_child_cohort_probabilities,
        MS_leaves_number_probabilities,
        dynT_user,
        dimT_user,
        GL_number,
        plantgen_core.DataCompleteness.SHORT,
        plantgen_core.DataCompleteness.SHORT,
        TT_hs_break,
        delais_TT_stop_del_axis,
        number_of_ears,
    )

    expected_phenT_abs = pandas.read_csv(default_expected_results_dir / "phenT_abs.csv")
    test_table_filepath = default_results.joinpath("phenT_abs.csv")
    phenT_abs.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        phenT_abs.values,
        expected_phenT_abs.values,
        relative_tolerance,
        absolute_tolerance,
    )

    expected_phenT = pandas.read_csv(default_expected_results_dir / "phenT.csv")
    test_table_filepath = default_results.joinpath("phenT.csv")
    phenT_.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        phenT_.values, expected_phenT.values, relative_tolerance, absolute_tolerance
    )

    expected_HS_GL_SSI_T = pandas.read_csv(
        default_expected_results_dir / "HS_GL_SSI_T.csv"
    )
    test_table_filepath = default_results.joinpath("HS_GL_SSI_T.csv")
    HS_GL_SSI_T.to_csv(
        test_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
    )
    print("The results have been saved to %s" % test_table_filepath)
    np.testing.assert_allclose(
        HS_GL_SSI_T.values,
        expected_HS_GL_SSI_T.values,
        relative_tolerance,
        absolute_tolerance,
    )


def test_gen_adel_input_data_from_min_min():
    reinit_random_state()
    dynT_user = pandas.read_csv(min_min_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(min_min_expected_results_dir / "dimT_user.csv")

    results = plantgen_interface.gen_adel_input_data(
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
    )

    expected_axeT = pandas.read_csv(min_min_expected_results_dir / "axeT.csv")
    expected_phenT_abs = pandas.read_csv(min_min_expected_results_dir / "phenT_abs.csv")
    expected_phenT = pandas.read_csv(min_min_expected_results_dir / "phenT.csv")
    expected_dynT = pandas.read_csv(min_min_expected_results_dir / "dynT.csv")
    expected_phenT_first = pandas.read_csv(
        min_min_expected_results_dir / "phenT_first.csv"
    )
    expected_HS_GL_SSI_T = pandas.read_csv(
        min_min_expected_results_dir / "HS_GL_SSI_T.csv"
    )
    expected_dimT = pandas.read_csv(min_min_expected_results_dir / "dimT.csv")
    expected_tilleringT = pandas.read_csv(
        min_min_expected_results_dir / "tilleringT.csv"
    )
    expected_cardinalityT = pandas.read_csv(
        min_min_expected_results_dir / "cardinalityT.csv"
    )

    to_compare = {
        "axeT": (expected_axeT, results[0]),
        "dimT": (expected_dimT, results[1]),
        "phenT": (expected_phenT, results[2]),
        "phenT_abs": (expected_phenT_abs, results[3]),
        "dynT": (expected_dynT, results[4]),
        "phenT_first": (expected_phenT_first, results[5]),
        "HS_GL_SSI_T": (expected_HS_GL_SSI_T, results[6]),
        "tilleringT": (expected_tilleringT, results[7]),
        "cardinalityT": (expected_cardinalityT, results[8]),
    }

    _check_results(
        to_compare,
        plantgen_core.DataCompleteness.MIN,
        plantgen_core.DataCompleteness.MIN,
    )


def test_gen_adel_input_data_from_short_short():
    reinit_random_state()
    dynT_user = pandas.read_csv(short_short_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(short_short_expected_results_dir / "dimT_user.csv")
    results = plantgen_interface.gen_adel_input_data(
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
    )

    expected_axeT = pandas.read_csv(short_short_expected_results_dir / "axeT.csv")
    expected_phenT_abs = pandas.read_csv(
        short_short_expected_results_dir / "phenT_abs.csv"
    )
    expected_phenT = pandas.read_csv(short_short_expected_results_dir / "phenT.csv")
    expected_dynT = pandas.read_csv(short_short_expected_results_dir / "dynT.csv")
    expected_phenT_first = pandas.read_csv(
        short_short_expected_results_dir / "phenT_first.csv"
    )
    expected_HS_GL_SSI_T = pandas.read_csv(
        short_short_expected_results_dir / "HS_GL_SSI_T.csv"
    )
    expected_dimT = pandas.read_csv(short_short_expected_results_dir / "dimT.csv")
    expected_tilleringT = pandas.read_csv(
        short_short_expected_results_dir / "tilleringT.csv"
    )
    expected_cardinalityT = pandas.read_csv(
        short_short_expected_results_dir / "cardinalityT.csv"
    )

    to_compare = {
        "axeT": (expected_axeT, results[0]),
        "dimT": (expected_dimT, results[1]),
        "phenT": (expected_phenT, results[2]),
        "phenT_abs": (expected_phenT_abs, results[3]),
        "dynT": (expected_dynT, results[4]),
        "phenT_first": (expected_phenT_first, results[5]),
        "HS_GL_SSI_T": (expected_HS_GL_SSI_T, results[6]),
        "tilleringT": (expected_tilleringT, results[7]),
        "cardinalityT": (expected_cardinalityT, results[8]),
    }

    _check_results(
        to_compare,
        plantgen_core.DataCompleteness.SHORT,
        plantgen_core.DataCompleteness.SHORT,
    )


def test_gen_adel_input_data_from_full_full():
    reinit_random_state()
    dynT_user = pandas.read_csv(full_full_expected_results_dir / "dynT_user.csv")
    dimT_user = pandas.read_csv(full_full_expected_results_dir / "dimT_user.csv")
    results = plantgen_interface.gen_adel_input_data(
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
    )

    expected_axeT = pandas.read_csv(full_full_expected_results_dir / "axeT.csv")
    expected_phenT_abs = pandas.read_csv(
        full_full_expected_results_dir / "phenT_abs.csv"
    )
    expected_phenT = pandas.read_csv(full_full_expected_results_dir / "phenT.csv")
    expected_dynT = pandas.read_csv(full_full_expected_results_dir / "dynT.csv")
    expected_phenT_first = pandas.read_csv(
        full_full_expected_results_dir / "phenT_first.csv"
    )
    expected_HS_GL_SSI_T = pandas.read_csv(
        full_full_expected_results_dir / "HS_GL_SSI_T.csv"
    )
    expected_dimT = pandas.read_csv(full_full_expected_results_dir / "dimT.csv")
    expected_tilleringT = pandas.read_csv(
        full_full_expected_results_dir / "tilleringT.csv"
    )
    expected_cardinalityT = pandas.read_csv(
        full_full_expected_results_dir / "cardinalityT.csv"
    )

    to_compare = {
        "axeT": (expected_axeT, results[0]),
        "dimT": (expected_dimT, results[1]),
        "phenT": (expected_phenT, results[2]),
        "phenT_abs": (expected_phenT_abs, results[3]),
        "dynT": (expected_dynT, results[4]),
        "phenT_first": (expected_phenT_first, results[5]),
        "HS_GL_SSI_T": (expected_HS_GL_SSI_T, results[6]),
        "tilleringT": (expected_tilleringT, results[7]),
        "cardinalityT": (expected_cardinalityT, results[8]),
    }

    _check_results(
        to_compare,
        plantgen_core.DataCompleteness.FULL,
        plantgen_core.DataCompleteness.FULL,
    )


def _check_results(to_compare, dynT_user_completeness, dimT_user_completeness):
    result_table_dir = tmp_results_directory.joinpath(
        "%s_%s" % (dynT_user_completeness.lower(), dimT_user_completeness.lower())
    )
    if not result_table_dir.exists():
        result_table_dir.mkdir()
    for key, value in to_compare.items():
        expected_table = value[0]
        result_table = value[1]
        result_table_filepath = result_table_dir.joinpath(key + ".csv")
        result_table.to_csv(
            result_table_filepath, na_rep="NA", index=False, float_format=FLOAT_FORMAT
        )
        print("The results have been saved to %s" % result_table_filepath)
        if "id_axis" in expected_table:
            assert (result_table["id_axis"] == expected_table["id_axis"]).all()
            result_table = result_table.drop(["id_axis"], axis=1)
            expected_table = expected_table.drop(["id_axis"], axis=1)
        np.testing.assert_allclose(
            result_table.values.astype(float),
            expected_table.values,
            relative_tolerance,
            absolute_tolerance,
        )


# Deactivate test as it modifies data : TO DO : pass tmp dir to node before executing
# def test_plantgen():
#     # test the visualea node
#     pm = PackageManager()
#     pm.init(verbose=False)
#     res = run(('openalea.adel.Tutorials', 'plantgen'), {}, pm=pm)
#     assert res == []


# if __name__ == '__main__':
#     test_init_axes()
#     test_phenology_functions()
#     test_plants_structure()
#     test_organs_dimensions()
#     test_axes_phenology()
#     test_gen_adel_input_data_from_min_min()
#     test_gen_adel_input_data_from_short_short()
#     test_gen_adel_input_data_from_full_full()
#     # test_plantgen()
