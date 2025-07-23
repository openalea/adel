from openalea.adel.dresser import (
    blade_dimension,
    stem_dimension,
    ear_dimension,
    dimension_table,
    AdelDress,
    AdelDressDyn,
)
from openalea.adel.data_samples import leaves


def test_dimension():
    blades = blade_dimension()
    stem = stem_dimension()
    ear = ear_dimension()
    dim = dimension_table()


def test_dresser():
    adel = AdelDress()
    g = adel.canopy()
    s = adel.scene(g)
    assert len(s) == 6
    stats = adel.midrib_statistics(g)
    assert all(stats["insertion_height"].values == (40, 50, 60))


def test_dresser_dyn():
    adel = AdelDressDyn()
    g = adel.canopy()


def test_mixture():
    l = leaves()
    adel = AdelDress(nplants=10, leaves={0: l[0], 1: l[0]}, species={0: 0.1, 1: 0.9})
    df = adel.canopy_table(adel.plant_references, adel.plant_species)
    assert 0.91 > df["species"].sum() * 1.0 / df["species"].size > 0.89
    g = adel.canopy()
    spec = [v for k, v in g.property("species").items() if k in g.vertices(1)]
    assert sum(spec) == 9
