import tempfile

import numpy as np
import pandas as pd
from pathlib import Path as path

from openalea.adel import postprocessing as pp

DATA_DIRPATH =  path(__file__).parent.resolve() / path("data/test_postprocessing")

INPUTS_DIRPATH = DATA_DIRPATH / "inputs"

OUTPUTS_DIRPATH = DATA_DIRPATH / "outputs"

ADEL_OUTPUT_FILENAME = "adel_output.csv"

RELATIVE_TOLERANCE = 10e-3
ABSOLUTE_TOLERANCE = 10e-3


# TO DO restore full test by validating new expectations
def test_aggregate_adel_output():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH / ADEL_OUTPUT_FILENAME)
    adel_output_df["species"] = "0"
    adel_output_aggregated_df = pp.aggregate_adel_output(
        adel_output_df, by=["TT", "species", "plant", "axe_id"]
    )
    adel_output_df.drop(["species"], axis=1, inplace=True)
    adel_output_aggregated_df.to_csv(
        OUTPUTS_DIRPATH / "actual_adel_aggregated_output.csv", index=False, na_rep="NA"
    )
    desired_adel_output_aggregated_df = pd.read_csv(
        OUTPUTS_DIRPATH / "desired_adel_aggregated_output.csv"
    )
    adel_output_aggregated_df = adel_output_aggregated_df.select_dtypes(
        include=[np.number]
    )
    desired_adel_output_aggregated_df = desired_adel_output_aggregated_df.select_dtypes(
        include=[np.number]
    )
    np.testing.assert_allclose(
        adel_output_aggregated_df.values,
        desired_adel_output_aggregated_df.values,
        RELATIVE_TOLERANCE,
        ABSOLUTE_TOLERANCE,
    )


def test_phenology():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH / ADEL_OUTPUT_FILENAME)
    adel_output_df["species"] = "0"
    phenology_df = pp.phenology(adel_output_df)
    phenology_df.drop(["species"], axis=1, inplace=True)
    phenology_df.to_csv(
        OUTPUTS_DIRPATH / "actual_phenology.csv", index=False, na_rep="NA"
    )
    desired_phenology_df = pd.read_csv(OUTPUTS_DIRPATH / "desired_phenology.csv")
    desired_phenology_df.drop(["has_ear"], axis=1, inplace=True)
    phenology_df = phenology_df.select_dtypes(include=[np.number])
    desired_phenology_df = desired_phenology_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(
        phenology_df.values,
        desired_phenology_df.values,
        RELATIVE_TOLERANCE,
        ABSOLUTE_TOLERANCE,
    )


def test_axis_statistics():
    adel_output_df = pd.read_csv(INPUTS_DIRPATH / ADEL_OUTPUT_FILENAME)
    adel_output_df["species"] = "0"
    axis_statistics_df, intermediate_df = pp.axis_statistics(
        adel_output_df, domain_area=1
    )
    axis_statistics_df.drop(["species"], axis=1, inplace=True)
    intermediate_df.drop(["species"], axis=1, inplace=True)
    axis_statistics_df.to_csv(
        OUTPUTS_DIRPATH / "actual_axis_statistics.csv", index=False, na_rep="NA"
    )
    intermediate_df.to_csv(
        OUTPUTS_DIRPATH / "actual_intermediate.csv", index=False, na_rep="NA"
    )

    desired_axis_statistics_df = pd.read_csv(
        OUTPUTS_DIRPATH / "desired_axis_statistics.csv"
    )
    desired_axis_statistics_df.drop(["has_ear"], axis=1, inplace=True)
    axis_statistics_df = axis_statistics_df.select_dtypes(include=[np.number])
    desired_axis_statistics_df = desired_axis_statistics_df.select_dtypes(
        include=[np.number]
    )
    np.testing.assert_allclose(
        axis_statistics_df.values,
        desired_axis_statistics_df.values,
        RELATIVE_TOLERANCE,
        ABSOLUTE_TOLERANCE,
    )

    desired_intermediate_df = pd.read_csv(OUTPUTS_DIRPATH / "desired_intermediate.csv")
    desired_intermediate_df.drop(["has_ear"], axis=1, inplace=True)
    intermediate_df = intermediate_df.select_dtypes(include=[np.number])
    desired_intermediate_df = desired_intermediate_df.select_dtypes(include=[np.number])
    np.testing.assert_allclose(
        intermediate_df.values,
        desired_intermediate_df.values,
        RELATIVE_TOLERANCE,
        ABSOLUTE_TOLERANCE,
    )


def test_plot_statistics():
    axis_statistics_df = pd.read_csv(INPUTS_DIRPATH / "axis_statistics.csv")
    axis_statistics_df["species"] = "0"
    plot_statistics_df = pp.plot_statistics(
        axis_statistics_df, plant_number=9, domain_area=1
    )
    plot_statistics_df.drop(["species"], axis=1, inplace=True)
    plot_statistics_df.to_csv(
        OUTPUTS_DIRPATH / "actual_plot_statistics.csv", index=False, na_rep="NA"
    )
    desired_plot_statistics_df = pd.read_csv(
        OUTPUTS_DIRPATH / "desired_plot_statistics.csv"
    )
    plot_statistics_df = plot_statistics_df.select_dtypes(include=[np.number])
    desired_plot_statistics_df = desired_plot_statistics_df.select_dtypes(
        include=[np.number]
    )
    np.testing.assert_allclose(
        plot_statistics_df.values,
        desired_plot_statistics_df.values,
        RELATIVE_TOLERANCE,
        ABSOLUTE_TOLERANCE,
    )


# if __name__ == '__main__':
#     test_aggregate_adel_output()
#     test_phenology()
#     test_axis_statistics()
#     test_plot_statistics()
