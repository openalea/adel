from openalea.adel.Stand import AgronomicStand
from openalea.adel.adelwheat_dynamic import AdelWheatDyn
from openalea.adel.AdelR import devCsv
from openalea.adel.plantgen_extensions import TillerEmission, TillerRegression, \
    AxePop, PlantGen, HaunStage
from openalea.adel.echap_leaf import echap_leaves
from openalea.mtg.mtg import MTG


def test_init_MTG_with_tillers():

    nplants=1
    sowing_density=250.
    plant_density=250.
    inter_row=0.15
    nff=12
    nsect=1
    seed=1
    leaves=echap_leaves(xy_model='Soissons_byleafclass')

    stand = AgronomicStand(sowing_density=sowing_density,
                           plant_density=plant_density, inter_row=inter_row,
                           noise=0.04, density_curve_data=None)
    em = TillerEmission(
        primary_tiller_probabilities={'T1': 1., 'T2': 1, 'T3': 1,
                                      'T4': 1, 'T5': 1, 'T6': 1, 'T7': 1})
    reg = TillerRegression(ears_per_plant=3)
    axp = AxePop(MS_leaves_number_probabilities={str(nff): 1}, Emission=em, Regression=reg)
    plants = axp.plant_list(nplants=nplants)
    hs = HaunStage(mean_nff=nff)
    pgen = PlantGen(HSfit=hs)
    axeT, dimT, phenT = pgen.adelT(plants)
    axeT = axeT.sort_values(['id_plt', 'id_cohort', 'N_phytomer'])
    devT = devCsv(axeT, dimT, phenT)
    adel = AdelWheatDyn(nplants=nplants, nsect=nsect, devT=devT, stand=stand,
                     seed=seed, sample='sequence', leaves=leaves,  scene_unit='m', devT_unit='cm')
    age = hs.TT(reg.hs_debreg(nff=nff))
    g = adel.setup_canopy(age)
    assert isinstance(adel, AdelWheatDyn)
    assert isinstance(g, MTG)
