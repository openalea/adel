"""Tutorial Reconstructing canopy from digitised data"""

import numpy

from openalea.adel.dresser import (
    blade_dimension,
    stem_dimension,
    ear_dimension,
    dimension_table,
    AdelDress,
)
from openalea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand

# input camille

blades = blade_dimension(
    length=[0] * 3 + [17.07, 12.1, 8.87],
    width=[0] * 3 + [1, 1.23, 1.01],
    ntop=[6, 5, 4, 3, 2, 1],
)
sheath = [0] * 3 + [11.21, 12.55, 13.59]
en = [0.55, 4.97, 8.28, 13.06, 13.36, 18.32]
h_ins = numpy.array(en).cumsum() + numpy.array(sheath)
stem = stem_dimension(
    ntop=[6, 5, 4, 3, 2, 1], h_ins=h_ins, d_internode=[0] * 3 + [0.34, 0.44, 0.31]
)
ear = ear_dimension(
    peduncle=14.8, ear=7.2, spike=8.8, projected_area_ear=10.28, d_peduncle=0.2
)
dimT = dimension_table(blades, stem, ear)

adel = AdelDress(dimT=dimT)
g = adel.canopy(nplants=1)
adel.plot(g)


# input Romain
#


blades = blade_dimension(
    length=[18.2, 21.1, 22.7, 17.4], area=[16, 22.8, 34, 34.6], ntop=[4, 3, 2, 1]
)
stem = stem_dimension(
    ntop=[4, 3, 2, 1],
    sheath=[11, 12.5, 14, 14.5],
    d_sheath=[0.2, 0.3, 0.4, 0.4],
    internode=[5, 8.6, 12.8, 18.6],
    d_internode=[0.2, 0.3, 0.3, 0.3],
)
ear = ear_dimension(peduncle=21.9, ear=9, projected_area_ear=15, d_peduncle=0.3)
dimT = dimension_table(blades, stem, ear)
# leaf shape database
leaves = Leaves()
stand = AgronomicStand(
    sowing_density=500, plant_density=500, inter_row=0.15, noise=0.03
)
adel = AdelDress(dimT=dimT, leaves=leaves, stand=stand)
g = adel.canopy(nplants=50)
adel.plot(g)


# mixture

leaves = {"soissons": Leaves(), "cap_horn": Leaves()}
stand = AgronomicStand(
    sowing_density=500, plant_density=500, inter_row=0.15, noise=0.03
)
adel = AdelDress(dimT=dimT, leaves=leaves, stand=stand)
g = adel.canopy(
    nplants=10,
    species={"soissons": 0.5, "cap_horn": 0.5},
    relative_inclination={"soissons": 1.5, "cap_horn": 1},
)
print([v for k, v in g.property("species").items() if k in g.vertices(1)])
print(adel.plot_statistics(g))
