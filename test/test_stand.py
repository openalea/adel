import tempfile

import numpy as np
import pandas

# from pathlib import Path as path
from pathlib import Path as path
from openalea.adel.stand.stand import post_processing

relative_tolerance = 10e-3
absolute_tolerance = 10e-3


def test_post_processing():
    adel_output_path = path(__file__).parent.resolve() / path("data/test_stand/adel_output.csv")

    global_postprocessing_file_path = path(__file__).parent.resolve() / path(
        tempfile.mktemp(prefix="global_postprocessing", suffix=".csv")
    )
    peraxis_postprocessing_file_path = path(__file__).parent.resolve() / path(
        tempfile.mktemp(prefix="peraxis_postprocessing", suffix=".csv")
    )

    # launch the postprocessing
    (
        global_postprocessing_file_path,
        peraxis_postprocessing_file_path,
        intermediate_results_file_path,
    ) = post_processing(
        adel_output_path,
        3,
        0.013824884792626729,
        global_postprocessing_file_path,
        peraxis_postprocessing_file_path,
    )

    expected_results_dir_path =  path(__file__).parent.resolve() / path("data/test_stand")
    # check intermediate results
    intermediate_results_array = pandas.read_csv(intermediate_results_file_path).values
    intermediate_results_array = np.delete(
        intermediate_results_array, (1, 3), 1
    ).astype(float)
    expected_intermediate_results_array = pandas.read_csv(
        expected_results_dir_path / "intermediate.csv"
    ).values
    expected_intermediate_results_array = np.delete(
        expected_intermediate_results_array, 2, 1
    ).astype(float)
    np.testing.assert_allclose(
        intermediate_results_array,
        expected_intermediate_results_array,
        relative_tolerance,
        absolute_tolerance,
    )

    # check per axis post processing results
    peraxis_postprocessing_results_array = pandas.read_csv(
        peraxis_postprocessing_file_path
    ).values
    peraxis_postprocessing_results_array = np.delete(
        peraxis_postprocessing_results_array, [0, 2], 1
    ).astype(float)
    expected_peraxis_postprocessing_results_array = pandas.read_csv(
        expected_results_dir_path / "peraxis_postprocessing.csv"
    ).values
    expected_peraxis_postprocessing_results_array = np.delete(
        expected_peraxis_postprocessing_results_array, [0, 2], 1
    ).astype(float)
    np.testing.assert_allclose(
        peraxis_postprocessing_results_array,
        expected_peraxis_postprocessing_results_array,
        relative_tolerance,
        absolute_tolerance,
    )

    # check global post-processing results
    global_postprocessing_results_array = pandas.read_csv(
        global_postprocessing_file_path
    ).values
    global_postprocessing_results_array = np.delete(
        global_postprocessing_results_array, 0, 1
    ).astype(float)
    expected_global_postprocessing_results_array = pandas.read_csv(
        expected_results_dir_path / "global_postprocessing.csv"
    ).values
    expected_global_postprocessing_results_array = np.delete(
        expected_global_postprocessing_results_array, 0, 1
    ).astype(float)
    np.testing.assert_allclose(
        global_postprocessing_results_array,
        expected_global_postprocessing_results_array,
        relative_tolerance,
        absolute_tolerance,
    )


# if __name__ == '__main__':
#     test_post_processing()
