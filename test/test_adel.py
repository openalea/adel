import pytest

import openalea.adel.data_samples as test_data
from openalea.adel.adel import Adel
import os


@pytest.fixture
def adel():
    return Adel()


@pytest.fixture()
def g():
    return test_data.adel_two_metamers()


def test_instantiate(adel):
    assert len(adel.leaves[0].xydb) == 6
    assert adel.stand.plant_density == 250
    assert adel.convUnit == 0.01
    assert adel.nplants == 1
    assert adel.domain_area == 1.0 / adel.stand.plant_density
    assert adel.positions[0] == (0, 0, 0)


def test_stand(adel):
    assert adel.aspect == "smart"
    assert adel.domain_area == 0.004
    assert adel.nplants == 1
    adel.new_stand(nplants=4, aspect="line")
    assert adel.domain_area == 0.016
    assert adel.nplants == 4
    adel.new_stand(nplants=4, aspect="square")
    assert adel.domain_area == adel.nplants * 1.0 / adel.stand.plant_density
    assert adel.nplants == 5
    adel.new_stand(nplants=4, aspect="smart")
    assert adel.domain_area == 0.016
    assert adel.nplants == 4


def test_get_axis(adel, g):
    assert g.nb_scales() == 6
    ms = adel.get_axis(g)
    assert ms.nb_scales() == 4


def test_scene(adel, g):
    s = adel.scene(g)
    assert len(s) == 6


def test_statistics(adel, g):
    areas = adel.get_exposed_areas(g)
    assert "green_area" in areas
    assert "species" in areas
    assert round(areas["green_area"].values[0], 2) == 0.31
    species = g.property("species")
    g.remove_property("species")
    areas = adel.get_exposed_areas(g)
    assert "species" in areas
    g.add_property("species")
    g.property("species").update(species)
    axstats = adel.axis_statistics(g)
    print(axstats["LAI totale"].round(2).values[0])
    assert axstats["LAI totale"].round(2).values[0] == 0.02
    pstats = adel.plot_statistics(g, axstats)
    assert pstats["Nbr.axe.tot.m2"][0] == 500


def test_midribs(adel, g):
    midribs = adel.get_midribs(g)
    midstats = adel.midrib_statistics(g)
    assert midstats["insertion_height"].values[0] == 2


def test_save_and_load(adel, g):
    fgeom, fg = adel.save(g)
    gg = adel.load()
    assert len(gg) == len(g)
    if os.path.exists(fgeom):
        os.remove(fgeom)
    if os.path.exists(fg):
        os.remove(fg)


def test_duplicated(adel, g):
    try:
        gg = adel.duplicated(g)
        assert False
    except ValueError:
        assert True
    adel.new_stand(nplants=2, duplicate=2)
    gg = adel.duplicated(g)
    assert gg.nb_vertices() == 1 + 2 * (g.nb_vertices() - 1)


def test_build_mtg(adel, g):
    pars = test_data.canopy_two_metamers()
    gg = adel.build_mtg(pars, stand=None)
    assert gg.nb_vertices() == g.nb_vertices()
