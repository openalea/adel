from openalea.adel.astk_interface import AdelWheat


def test_static():
    nplants = 1
    adel = AdelWheat(nplants=nplants)
    g = adel.setup_canopy(age=100)
    assert len(g.vertices()) > 20
    assert list(g.property("geometry").values())[0].isValid()


def test_statistics():
    nplants = 1
    adel = AdelWheat(nplants=nplants)
    g = adel.setup_canopy(age=100)
    areas = adel.get_exposed_areas(g)
    assert "green_area" in areas
    pstats = adel.plot_statistics(g)


# to be repaired
# def test_dynamic():
#     nplants = 1
#     adel = AdelWheat(nplants=nplants)
#     g = adel.setup_canopy(age=100)
#     timing = [TimeControlSet(dt=100) for _ in range(2)]
#     for tc in timing:
#         adel.grow(g,tc)
#     return g
