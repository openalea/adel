from openalea.adel.astk_interface import AdelWheat
from openalea.astk.Weather import sample_weather


seq, weather = sample_weather()
wdata = weather.get_weather(seq)

adel = AdelWheat(nsect=2)

g = adel.setup_canopy(100)
adel.grow(g, wdata)
