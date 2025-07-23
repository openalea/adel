from itertools import cycle
from math import pi, sqrt
from operator import itemgetter
from pathlib import Path as path
from random import random, sample, normalvariate

import numpy as np
from numpy import cos, linspace, sin
from numpy.random import vonmises
import pandas as pd

import openalea.adel.postprocessing as pp


def regular(nb_plants, nb_rank, dx, dy, nx=None):
    if nx is None:
        nx = int(nb_plants / nb_rank)
    ny = nb_rank
    domain = ((0, 0), (nx * dx, ny * dy))
    return [
        (i * dx + dx / 2.0, j * dy + dy / 2.0, 0.0)
        for j in range(ny)
        for i in range(nx)
    ], domain


def randomise_position(position, radius):
    az = random() * 2 * np.pi
    r = random() * radius
    dx = r * cos(az)
    dy = r * sin(az)
    x, y, z = position
    return x + dx, y + dy, z


def regular_plot(
    inter_plant,
    inter_row,
    nrow,
    plant_per_row,
    noise=0,
    convunit=100,
    center_scene=True,
):
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = regular(
        nrow * plant_per_row, int(nrow), dx, dy, int(plant_per_row)
    )
    domain_area = (
        abs(domain[1][0] - domain[0][0])
        / convunit
        * abs(domain[1][1] - domain[0][1])
        / convunit
    )

    # sorting by ranks
    positions = sorted(positions, key=itemgetter(1, 0))
    # add noise
    if noise > 0:
        positions = [randomise_position(x, noise * convunit) for x in positions]
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x, y, z in positions]
        domain = (
            (domain[0][0] - xc, domain[0][1] - yc),
            (domain[1][0] - xc, domain[1][1] - yc),
        )

    return positions, domain, domain_area


def agronomicplot(
    length,
    width,
    sowing_density,
    plant_density,
    inter_row,
    noise=0,
    convunit=100,
    center_scene=True,
):
    """Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
    length (m) is plot dimension along row direction
    width (m) is plot dimension perpendicular to row direction
    sowing density is the density of seeds sawn
    plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
    inter_row (m) is for the  distance between rows
    noise (m) is the radius of the circle where individual plant are randomly positionned
    convunit is the conversion factor from meter to scene unit
    center_scene allows to center the position arround origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain

    Rows are parrallel to x-axis
    Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
    """

    # find a (nrow, plant_per_row) sowing design that best match plot dimension
    inter_plant = 1.0 / inter_row / sowing_density
    nrow = max(1, int(round(float(width) / inter_row)))
    plant_per_row = max(1, int(round(float(length) / inter_plant)))
    positions, domain, domain_area = regular_plot(
        inter_plant,
        inter_row,
        nrow,
        plant_per_row,
        noise=noise,
        convunit=convunit,
        center_scene=center_scene,
    )
    n_emerged = int(round(len(positions) * plant_density / sowing_density))
    positions = sample(positions, n_emerged)

    return n_emerged, positions, domain, domain_area, 1.0 / convunit


def agronomicplotwithdistributions(
    length,
    width,
    sowing_density,
    plant_density,
    inter_row,
    mu=0.0,
    kappa=3.0,
    noise=0,
    convunit=100,
):
    """
    Returns the
        - number of plants
        - the positions
        - azimuthal orientation of each plant
        - the domain
        - the simulated density
        of a micro-plot specified with agronomical variables

    Inputs:
        - length (m) is plot dimension along row direction
        - width (m) is plot dimension perpendicular to row direction
        - sowing density is the density of seeds sawn
        - plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
        - inter_row (m) is for the  distance between rows
        - mu, kappa are the von Mises parameters of the azimuthal distribution
        - noise (%), indicates the precision of the sowing for the inter plant spacing
        - unit (m or cm) is for the unit of the position and domain
    """
    inter_plant = 1.0 / inter_row / sowing_density
    nrow = max(1, int(float(width) / inter_row))
    plant_per_row = max(1, int(float(length) / inter_plant))
    nplants = nrow * plant_per_row
    positions, domain = regular(
        nplants, nrow, inter_plant * convunit, inter_row * convunit
    )
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = sample(positions, n_emerged)
    density = int(
        n_emerged
        / (
            abs(domain[1][0] - domain[0][0])
            / convunit
            * abs(domain[1][1] - domain[0][1])
            / convunit
        )
    )
    azimuths = vonmises(mu, kappa, nplants)
    return n_emerged, positions, azimuths, domain, density


def regularband(nb_plants, nb_rank, dx, dy):
    nx = int(nb_plants / nb_rank)
    ny = nb_rank
    db = sample(range(100), nb_plants)
    domain = ((0, 0), (nx * dx, ny * dy))
    return [
        (
            i * dx + dx / 2.0,
            j * dy + dy / 2.0 + (db[i + j] - 50) / 100.0 * dy / 3.0,
            0.0,
        )
        for j in range(nb_rank)
        for i in range(nx)
    ], domain


def concentric(nb_plants, distance_plant):
    nb_circle = int(sqrt(nb_plants))
    dp = distance_plant
    points = []
    for i in range(nb_circle):
        n = 2 * i + 1
        theta = linspace(0, pi, n)
        radius = i * dp
        x = cos(theta) * radius
        y = sin(theta) * radius
        points.extend([(x[i], y[i], 0.0) for i in range(len(x))])

    return points


def uniform(density, nb_plants):
    return []


def planter(scene, distribution):
    """
    Returns a CanestraScene by positionning all elements.
    """
    from openalea.plantgl.all import AxisRotation, Transform4, Translation, Vector3

    from openalea.adel import symbol

    new_scene = symbol.CanScene()

    def transfo(plant, pt):
        r = AxisRotation((0, 0, 1), random() * pi).getMatrix()
        return plant.transform2(Transform4(Translation(pt).getMatrix() * r))

    l = [
        new_scene.add_plant(p)
        for p in map(transfo, cycle(iter(scene.values())), iter(distribution))
    ]
    return new_scene


def sample_selection(points, gap_fraction):
    """
    Choose a sample from a list of points.
    The gap fraction indicates the proportion of points to remove.
    """
    if gap_fraction == 0:
        return points
    n = len(points)
    k = int(n * (1.0 - gap_fraction / 100.0))
    return sample(points, k)


def clumping_selection(points, nb_clumps, radius_mu, radius_sigma):
    """
    Select a set of points based on andomly choosen groups
    with a radius following a normal distribution.
    """
    n = len(points)
    centers = sample(list(range(n)), nb_clumps)

    for _ in centers:
        radius = normalvariate(radius_mu, radius_sigma)
        r2 = radius**2
        # filter to select all the poit in a ball of radius r2
        # Add these points to the result


def sample_regular_gaps(points, pattern=[0, 1]):
    """Sample points along pattern.
    Pattern is replicated to become as long as points and then used as a filter (0= gap) on points
    Returns positions of plants and gaps
    """

    if not isinstance(points, list):
        points = [points]

    p = pattern
    length = len(points)
    p = p * (length / len(p)) + [p[i] for i in range(length % len(p))]

    # selection = compress(points,p) only python >= 2.7
    return [point for point, i in zip(points, p) if i], [
        point for point, i in zip(points, p) if not i
    ]


def post_processing(
    adel_output_path="",
    plant_number=0,
    domain_area=0.0,
    global_postprocessing_path="",
    peraxis_postprocessing_path="",
    convUnit=0.01,
):
    """
    Apply post processing on the ADEL output.

    For one TT, the variables calculated are : TODO: list and explain the
    variables.
    The results of the post processing are added at the end of two csv files: one
    for the global results, and another one for the results per axis. These files
    are created if they don't exist yet.
    Furthermore, for each call to the current function, intermediate post
    processing results are saved in a new temporary file.

    :Parameters:

        - `adel_output_path` (str) -
          path to the csv file describing ADEL output. This file contains data
          for one TT.

        - `plant_number` (int) -
          the number of plants simulated by ADEL

        - `domain_area` (float) -
          the area of the domain on which ADEL simulation has been done.

        - `global_postprocessing_path` (str) -
          path to the csv file describing the results of the post processing
          for each TT.

        - `peraxis_postprocessing_path` (str) -
          path to the csv file describing the results, for each TT, of the post
          processing for each axis

        - `convUnit` (float) -
          Factor to convert the length unit of adel output to meter. Default '0.01'.

    :Returns:
        Three paths:
            * one for the global post processing results,
            * one for the post processing results per axis,
            * and one for the intermediate results.

    :Returns Type:
        tuple of 3 str

    .. warning:: adel_output_path, global_postprocessing_path and
                 peraxis_postprocessing_path must be non-empty,
                 plant_number and domain_area must be non-null.

    """

    if adel_output_path == "":
        raise Exception("adel_output_path == ''")
    if plant_number == 0:
        raise Exception("plant_number == 0")
    if domain_area == 0.0:
        raise Exception("domain_area == 0.0")
    if global_postprocessing_path == "":
        raise Exception("global_postprocessing_path == ''")
    if peraxis_postprocessing_path == "":
        raise Exception("peraxis_postprocessing_path == ''")
    if convUnit == 0.0:
        raise Exception("convUnit == 0.0")

    adel_output_path = path(adel_output_path)
    adel_output_df = pd.read_csv(adel_output_path)
    adel_output_df["species"] = "0"
    aggregated_adel_output_df = pp.aggregate_adel_output(
        adel_output_df, by=["plant", "axe_id"]
    )
    phenology_df = pp.phenology(adel_output_df)
    phenology_df.drop(["species"], axis=1, inplace=True)
    axis_statistics_df = pp.axis_statistics(adel_output_df, domain_area, convUnit)[0]
    plot_statistics_df = pp.plot_statistics(
        axis_statistics_df, plant_number, domain_area
    )
    axis_statistics_df.drop(["species"], axis=1, inplace=True)
    plot_statistics_df.drop(["species"], axis=1, inplace=True)
    intermediate_df = pd.merge(
        aggregated_adel_output_df, phenology_df, on=["TT", "plant", "axe_id"]
    )

    import tempfile

    TT = adel_output_df["TT"][0]  # there is only one TT
    intermediate_file_suffix = "-%d.csv" % int(TT)
    intermediate_path = path(tempfile.mktemp(intermediate_file_suffix))
    intermediate_df.to_csv(intermediate_path, na_rep="NA", index=False)

    axis_statistics_df.insert(0, "Filename", intermediate_path.name)
    peraxis_postprocessing_path_ = path(peraxis_postprocessing_path)
    if peraxis_postprocessing_path_.exists():
        old_peraxis_postprocessing_df = pd.read_csv(peraxis_postprocessing_path_)
        axis_statistics_df = pd.concat(
            [old_peraxis_postprocessing_df, axis_statistics_df], ignore_index=True
        )
    axis_statistics_df.to_csv(peraxis_postprocessing_path_, na_rep="NA", index=False)

    plot_statistics_df.insert(0, "Filename", intermediate_path.name)
    global_postprocessing_path_ = path(global_postprocessing_path)
    if global_postprocessing_path_.exists():
        old_global_postprocessing_df = pd.read_csv(global_postprocessing_path_)
        plot_statistics_df = pd.concat(
            [old_global_postprocessing_df, plot_statistics_df], ignore_index=True
        )
    plot_statistics_df.to_csv(global_postprocessing_path_, na_rep="NA", index=False)

    return (
        global_postprocessing_path,
        peraxis_postprocessing_path,
        str(intermediate_path),
    )
