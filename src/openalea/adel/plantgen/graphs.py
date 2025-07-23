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
Routines to plot the outputs of :mod:`openalea.adel.plantgen`.

Authors: M. Abichou, B. Andrieu, C. Chambon
"""

import os

from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd

from openalea.adel.plantgen import params


def plot_HS_GL_SSI_T(
    HS_GL_SSI_T,
    dynT,
    config,
    id_phen_to_plot=None,
    dynamics_to_plot=None,
    plot_most_frequent_main_stem=True,
    plot_non_regressive_tillers=True,
    plot_regressive_tillers=True,
    plot_filepath=None,
):
    """
    Plot the HS, GL and SSI of `id_phen_to_plot`.

    :Parameters:

        - `HS_GL_SSI_T` (:class:`pd.DataFrame`) - the table HS_GL_SSI_T.
        - `dynT` (:class:`pd.DataFrame`) - the table dynT.
        - `config` (:class:`dict`) - a dictionary which stores the configuration used for the construction. `config` must contain at least the GL decimal
          numbers measured at several thermal times (including the senescence end).
        - `id_phen_to_plot` (:class:`list`) - the list of id_phen to plot. If None (the default) or empty, then plot all the id_phen.
        - `dynamics_to_plot` (:class:`list`) - the list of dynamics to plot. The available dynamic are: 'HS', 'GL' and 'SSI'.
          If None (the default) or empty, then plot all the dynamics.
        - `plot_most_frequent_main_stem` (:class:`bool`) - whether to plot the most frequent main stem or not. If true, dynT must not be `None`.
        - `plot_non_regressive_tillers` (:class:`bool`) - whether to plot the non regressive tillers or not. Non regressive tillers have id_dim ending by '1'. Default is to plot the non regressive tillers.
        - `plot_regressive_tillers` (:class:`bool`) - whether to plot the regressive tillers or not. Regressive tillers have id_dim ending by '0'. Default is to plot the regressive tillers.
        - `plot_filepath` (:class:`str`) - the path of the file to save the plot in.
          If `None` (the default), do not save the plot but display it.

    :Examples:

        # plot HS, GL and SSI of id_phen 1101, 1111 and 4081
        import pandas as pd
        HS_GL_SSI_T = pd.read_csv('HS_GL_SSI_T.csv')
        dynT = pd.read_csv('dynT.csv')
        config = {'GL_number': {1117.0: 5.6, 1212.1:5.4, 1368.7:4.9, 1686.8:2.4, 1880.0:0.0}}
        from openalea.adel.plantgen import graphs
        graphs.plot_HS_GL_SSI_T(HS_GL_SSI_T, dynT, config, id_phen_to_plot=[4081, 1111, 1101])

    """
    HS_GL_SSI_T_to_plot = HS_GL_SSI_T.copy()

    id_cohort, N_phytomer_potential, t0, t1, TT_flag_ligulation, n0, n1, n2 = dynT.loc[
        dynT.first_valid_index(),
        [
            "id_cohort",
            "N_phytomer_potential",
            "t0",
            "t1",
            "TT_flag_ligulation",
            "n0",
            "n1",
            "n2",
        ],
    ]
    most_frequent_main_stem_id_phen = int(
        "".join([str(int(id_cohort)), str(int(N_phytomer_potential)).zfill(2), "1"])
    )

    if id_phen_to_plot is None or len(id_phen_to_plot) == 0:
        id_phen_to_plot = HS_GL_SSI_T_to_plot.id_phen.unique()
    else:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T_to_plot[
            HS_GL_SSI_T_to_plot["id_phen"].isin(id_phen_to_plot)
        ]

    if dynamics_to_plot is None or len(dynamics_to_plot) == 0:
        dynamics_to_plot = ["HS", "GL", "SSI"]
    else:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T_to_plot[["id_phen", "TT"] + dynamics_to_plot]

    if not plot_non_regressive_tillers:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T_to_plot[
            HS_GL_SSI_T_to_plot.id_phen.astype(str).str[-1].astype(int) != 1
        ]

    if not plot_regressive_tillers:
        HS_GL_SSI_T_to_plot = HS_GL_SSI_T_to_plot[
            HS_GL_SSI_T_to_plot.id_phen.astype(str).str[-1].astype(int) != 0
        ]

    LINE_STYLES = ["-", "--", "-.", ":"]

    DYNAMIC_TO_COLOR_MAPPING = {"HS": "b", "GL": "g", "SSI": "r"}

    plt.figure()
    plot_ = plt.subplot(111)

    axis_num = 0
    for id_phen, group in HS_GL_SSI_T_to_plot.groupby("id_phen"):
        axis_num += 1
        line_style = LINE_STYLES[axis_num % len(LINE_STYLES)]
        for dynamic_to_plot in dynamics_to_plot:
            plot_.plot(
                group.TT,
                group[dynamic_to_plot],
                linestyle=line_style,
                color=DYNAMIC_TO_COLOR_MAPPING[dynamic_to_plot],
                label="{} - {}".format(id_phen, dynamic_to_plot),
            )

        if plot_most_frequent_main_stem and id_phen == most_frequent_main_stem_id_phen:
            plot_.plot(
                [t0, t1, TT_flag_ligulation] + list(config["GL_number"].keys()),
                [n0, n1, n2] + list(config["GL_number"].values()),
                linestyle="",
                marker="D",
                color="k",
                label="{} - {}".format(id_phen, "measured data"),
            )

    plot_.set_xlabel("Thermal time")
    plot_.set_ylabel("Decimal leaf number")
    plot_.legend(prop={"size": 10}, framealpha=0.5)

    plot_.set_title("{} - {}".format(tuple(id_phen_to_plot), tuple(dynamics_to_plot)))

    xmin, xmax = plot_.get_xlim()
    x_margin = (xmax - xmin) / 100.0
    plot_.set_xlim(xmin - x_margin, xmax + x_margin)
    ymin, ymax = plot_.get_ylim()
    y_margin = (ymax - ymin) / 100.0
    plot_.set_ylim(ymin - y_margin, ymax + y_margin)
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format="PNG")
        plt.close()


def plot_dimT(
    dimT,
    measured_id_dim=None,
    relative_index_phytomer=False,
    dimensions_to_plot=None,
    id_dim_to_plot=None,
    id_cohort_to_plot=None,
    plot_non_regressive_tillers=True,
    plot_regressive_tillers=True,
    plots_dirpath=None,
):
    """
    Plot the dimensions in `dimT` according to filters `dimensions_to_plot`, `id_dim_to_plot`, `id_cohort_to_plot`, `plot_non_regressive_tillers` and `plot_regressive_tillers`.

    :Parameters:

        - `dimT` (:class:`pd.DataFrame`) - the table dimT.
        - `measured_id_dim` (:class:`list`) - the list of id_dim for which we have measured data. For these id_dim we plot no line but thick black markers.
          None (the default) or empty means that we do not represent measured data.
        - `relative_index_phytomer` (:class:`bool`) - if True: display the index relative to the phytomers of the main stem.
          If False (the default), display the absolute index of the phytomers.
        - `dimensions_to_plot` (:class:`list`) - the list of dimensions to plot. If None (the default) or empty, then plot all the dimensions.
        - `id_dim_to_plot` (:class:`list`) - the list of id_dim to plot. If None (the default) or empty, then plot all the id_dim.
        - `id_cohort_to_plot` (:class:`list`) - the list of id_cohort to plot. If None (the default) or empty, then plot all the id_cohort.
        - `plot_non_regressive_tillers` (:class:`bool`) - whether to plot the non regressive tillers or not. Non regressive tillers have id_dim ending by '1'. Default is to plot the non regressive tillers.
        - `plot_regressive_tillers` (:class:`bool`) - whether to plot the regressive tillers or not. Regressive tillers have id_dim ending by '0'. Default is to plot the regressive tillers.
        - `plots_dirpath` (:class:`str`) - the path of the directory to save the plots in.
          If `None`, do not save the plots but display them.

    :Examples:

        # plot L_blade
        import pandas as pd
        dimT = pd.read_csv('dimT.csv')
        from openalea.adel.plantgen import graphs
        graphs.plot_dimT(dimT, dimensions_to_plot=['L_blade'])

    """

    DIM_T_KEY = ["id_dim", "index_phytomer"]

    MARKERS = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != " ":
                MARKERS.append(m)
        except TypeError:
            pass

    MARKERS.extend(
        [
            r"$\lambda$",
            r"$\bowtie$",
            r"$\circlearrowleft$",
            r"$\clubsuit$",
            r"$\checkmark$",
        ]
    )

    COLORS = ("b", "g", "r", "c", "m", "y")

    if dimensions_to_plot is None or len(dimensions_to_plot) == 0:
        dimT_to_plot = dimT
        dimensions_to_plot = dimT.columns.difference(DIM_T_KEY)
    else:
        dimT_to_plot = dimT[DIM_T_KEY + dimensions_to_plot]

    if id_dim_to_plot is None or len(id_dim_to_plot) == 0:
        dimT_to_plot = dimT_to_plot
    else:
        dimT_to_plot = dimT_to_plot[dimT_to_plot.id_dim.isin(id_dim_to_plot)]

    if id_cohort_to_plot is not None and len(id_cohort_to_plot) != 0:
        dimT_to_plot = dimT_to_plot[
            dimT_to_plot.id_dim.astype(str).str[0].astype(int).isin(id_cohort_to_plot)
        ]

    if not plot_non_regressive_tillers:
        dimT_to_plot = dimT_to_plot[
            dimT_to_plot.id_dim.astype(str).str[-1].astype(int) != 1
        ]

    if not plot_regressive_tillers:
        dimT_to_plot = dimT_to_plot[
            dimT_to_plot.id_dim.astype(str).str[-1].astype(int) != 0
        ]

    xlabel = "index_phytomer"
    if relative_index_phytomer:
        xlabel = "relative_" + xlabel

    for dimension in dimensions_to_plot:
        plt.figure()
        current_axis = plt.subplot(111)
        dimension_to_plot = dimT_to_plot[DIM_T_KEY + [dimension]]
        grouped = dimension_to_plot.groupby("id_dim")
        axis_num = 0
        lines_for_the_legend = []
        id_dims_for_the_legend = []
        for id_dim, group in grouped:
            axis_num += 1
            if relative_index_phytomer:
                id_cohort = int(str(int(id_dim))[:-3])
                if id_cohort == 1:  # MS
                    index_phytomer_to_plot = group.index_phytomer
                else:
                    index_phytomer_to_plot = group.index_phytomer - (
                        params.SLOPE_SHIFT_MS_TO_TILLERS * id_cohort
                    )
            else:
                index_phytomer_to_plot = group.index_phytomer
            if measured_id_dim is not None and id_dim in measured_id_dim:
                color = "k"
                markersize = 10
                linestyle = ""
            else:
                color = COLORS[axis_num % len(COLORS)]
                markersize = 7
                linestyle = "-"

            marker = MARKERS[axis_num % len(MARKERS)]
            dimension_to_plot = group[dimension]

            current_axis.plot(
                index_phytomer_to_plot,
                dimension_to_plot,
                color=color,
                linestyle="",
                marker=marker,
                markersize=markersize,
            )
            interpolated_index_phytomer_to_plot = np.linspace(
                index_phytomer_to_plot.min(), index_phytomer_to_plot.max(), num=1000
            )
            index_phytomer_to_plot_length = len(index_phytomer_to_plot)
            if index_phytomer_to_plot_length == 1:
                interpolated_dimension_to_plot = dimension_to_plot.repeat(
                    len(interpolated_index_phytomer_to_plot)
                )
            else:
                interpolated_dimension_to_plot = (
                    interpolate.InterpolatedUnivariateSpline(
                        index_phytomer_to_plot,
                        dimension_to_plot,
                        k=min(index_phytomer_to_plot_length - 1, 5),
                    )(interpolated_index_phytomer_to_plot)
                )
            current_axis.plot(
                interpolated_index_phytomer_to_plot,
                interpolated_dimension_to_plot,
                color=color,
                linestyle=linestyle,
            )
            line_for_the_legend = plt.Line2D(
                (0, 1), (0, 0), color=color, marker=marker, linestyle=linestyle
            )
            lines_for_the_legend.append(line_for_the_legend)
            id_dims_for_the_legend.append(str(int(id_dim)))

        current_axis.set_xlabel(xlabel)
        current_axis.set_ylabel(dimension)
        current_axis.legend(
            lines_for_the_legend,
            id_dims_for_the_legend,
            prop={"size": 10},
            framealpha=0.5,
        )
        current_axis.set_title(dimension)
        xmin, xmax = current_axis.get_xlim()
        x_margin = (xmax - xmin) / 100.0
        current_axis.set_xlim(xmin - x_margin, xmax + x_margin)
        ymin, ymax = current_axis.get_ylim()
        y_margin = (ymax - ymin) / 100.0
        current_axis.set_ylim(ymin - y_margin, ymax + y_margin)

        if plots_dirpath is None:
            plt.show()
        else:
            plt.savefig(
                os.path.join(plots_dirpath, "{}.png".format(dimension)),
                dpi=200,
                format="PNG",
            )
            plt.close()


def plot_tillering_dynamic(axeT, plants_density, TT_step=10, plots_dirpath=None):
    """
    Plot the dynamic of tillering, i.e. the evolution of the density of active axes when TT varies.

    A non regressive axis is active at TT if TT >= TT_em_phytomer1.
    A regressive axis is active at TT if TT_em_phytomer1 <= TT < TT_stop_axis.

    :Parameters:

        - `axeT` (:class:`pd.DataFrame`) - the table axeT.
        - `plants_density` (:class:`int`) - the number of plants per square meter.
        - `TT_step` (:class:`int`) - the thermal time step of the plot. Default is 10.
        - `plots_dirpath` (:class:`str`) - the path of the directory to save the plots in.
          If `None`, do not save the plots but display them.

    :Examples:

        # plot tillering dynamic
        import pandas as pd
        axeT = pd.read_csv('axeT.csv')
        from openalea.adel.plantgen import graphs
        graphs.plot_tillering_dynamic(axeT, plants_density=250)

    """

    number_of_plants = axeT.id_plt.nunique()
    domain_area = number_of_plants / float(plants_density)

    xmin = axeT.TT_em_phytomer1.min() - TT_step
    xmax = axeT.TT_stop_axis.max() + TT_step
    if TT_step >= 10:
        xmin = round(xmin, -1)
        xmax = round(xmax, -1)
    TT_grid = np.arange(xmin, xmax + TT_step, TT_step)

    axeT_non_regressive = axeT[axeT.id_ear.notnull()]
    null_id_ear_index = axeT.index.difference(axeT_non_regressive.index)
    axeT_regressive = axeT.loc[null_id_ear_index, :]

    densities_of_active_axes_per_square_meter = []
    densities_of_active_axes_per_plant = []

    for TT in TT_grid:
        number_of_non_regressive_active_axes = len(
            axeT_non_regressive[TT >= axeT_non_regressive.TT_em_phytomer1]
        )
        number_of_regressive_active_axes = len(
            axeT_regressive[
                (TT >= axeT_regressive.TT_em_phytomer1)
                & (TT < axeT_regressive.TT_stop_axis)
            ]
        )
        total_number_of_active_axes = (
            number_of_non_regressive_active_axes + number_of_regressive_active_axes
        )
        densities_of_active_axes_per_square_meter.append(
            total_number_of_active_axes / domain_area
        )
        densities_of_active_axes_per_plant.append(
            total_number_of_active_axes / float(number_of_plants)
        )

    for densities_of_active_axes, type_of_density, file_suffix in (
        (
            densities_of_active_axes_per_square_meter,
            "per square meter",
            "per_square_meter",
        ),
        (densities_of_active_axes_per_plant, "per plant", "per_plant"),
    ):
        plt.figure()
        current_axis = plt.subplot(111)
        densities_of_active_axes_series = pd.Series(
            densities_of_active_axes, index=TT_grid
        )
        current_axis = densities_of_active_axes_series.plot()
        current_axis.set_xlabel("Thermal time")
        current_axis.set_ylabel("Density of active axes " + type_of_density)
        current_axis.set_title("Tillering dynamic")
        xmin, xmax = current_axis.get_xlim()
        x_margin = (xmax - xmin) / 100.0
        current_axis.set_xlim(xmin - x_margin, xmax + x_margin)
        ymin, ymax = current_axis.get_ylim()
        y_margin = (ymax - ymin) / 100.0
        current_axis.set_ylim(ymin - y_margin, ymax + y_margin)
        if plots_dirpath is None:
            plt.show()
        else:
            plt.savefig(
                os.path.join(
                    plots_dirpath, "{}_{}.png".format("tillering_dynamic", file_suffix)
                ),
                dpi=200,
                format="PNG",
            )
            plt.close()


def plot_tiller_probabilities(dynT, config, plot_filepath=None):
    """
    Plot the probability of each tiller.

    :Parameters:

        - `dynT` (:class:`pd.DataFrame`) - the table dynT.
        - `config` (:class:`int`) - `config` (:class:`dict`) - a dictionary which stores the configuration used for the construction.
          `config` must contain at least the number of plants.
        - `plot_filepath` (:class:`str`) - the path of the file to save the plot in.
          If `None` (the default), do not save the plot but display it.

    :Examples:

        # plot tiller probabilities
        import pandas as pd
        dynT = pd.read_csv('dynT.csv')
        from openalea.adel.plantgen import graphs
        graphs.plot_tiller_probabilities(dynT, config)

    """
    summed_cardinalities = dynT[["id_axis", "cardinality"]].groupby("id_axis").sum()
    hist = summed_cardinalities.cardinality.values
    normalized_hist = hist.astype(float) / config["plants_number"]
    bin_edges = np.arange(1, len(summed_cardinalities) + 1)
    width = 0.5
    plt.bar(bin_edges, normalized_hist, width=width, align="center")
    plt.title("Axis probabilities")
    plt.xlabel("id_axis")
    plt.ylabel("Probability")
    ymin, ymax = plt.ylim()
    y_margin = (ymax - ymin) / 10.0
    plt.ylim(ymin, ymax + y_margin)
    plt.xticks(bin_edges, summed_cardinalities.index)
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format="PNG")
        plt.close()
