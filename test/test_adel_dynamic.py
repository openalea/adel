"""Test use cases of adel/fspm wheat coupling"""

from openalea.adel.adel_dynamic import AdelDyn
import openalea.adel.data_samples as test_data


def test_build_stand():
    adel = AdelDyn()
    axeT = test_data.axeTable()
    g = adel.build_stand(axeT)
    assert g.nb_vertices(scale=1) == 1
    assert g[2]["hasEar"]


def test_add_metamer():
    adel = AdelDyn()
    axeT = test_data.axeTable()
    phytoT = test_data.phytoT()
    g = adel.build_stand(axeT)
    vid = adel.add_metamer(g, phytoT)
    labels = g.property("label")
    new_metamer = g.node(vid)
    metamers = [vid for vid in labels if labels[vid].startswith("metamer")]
    assert len(metamers) == 2
    internode, sheath, blade = new_metamer.components()
    assert blade.length == 0
    elts = [c.label for c in blade.components()]
    assert elts == ["baseElement", "topElement"]


def test_build_mtg():
    pars = test_data.canopy_two_metamers()
    adel = AdelDyn()
    g = adel.build_mtg(pars, stand=None)


def test_update_geometry():
    adel = AdelDyn()
    axeT = test_data.axeTable()
    phytoT = test_data.phytoT()
    g = adel.build_stand(axeT)
    vid = adel.add_metamer(g, phytoT)
    new_metamer = g.node(vid)
    internode, sheath, blade = new_metamer.components()
    assert blade.length == 0
    # grow leaf and check components
    blade.length = 6
    blade.visible_length = 3
    adel.update_geometry(g)
    assert blade.area > 0
    elts = [c.label for c in blade.components()]
    assert len(elts) > 2
