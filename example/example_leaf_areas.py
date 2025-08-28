"""Check progression of leaf areas"""

import matplotlib.pyplot as plt

plt.ion()

# Imports for wheat
try:
    import pickle as pickle
except:
    import pickle
from openalea.adel.newmtg import move_properties, adel_ids
from openalea.echap.architectural_reconstructions import EchapReconstructions
from openalea.alep.disease_outputs import initiate_all_adel_septo_recorders

# Imports for weather
from alinea.astk.TimeControl import *
from alinea.echap.weather_data import *
from alinea.alep.alep_weather import linear_degree_days

# Temporary
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer

from random import seed

seed(0)

# Initialize wheat plant
nsect = 5
reconst = EchapReconstructions()
adel = reconst.get_reconstruction(
    name="Tremie12", nplants=1, nsect=nsect, disc_level=5, aspect="line"
)
g = adel.setup_canopy(age=300.0)

# Manage weather and time control
start_date = "2012-03-01 12:00:00"
end_date = "2012-08-01 01:00:00"
weather = Boigneville_2011_2012()
weather.check(
    varnames=["degree_days"],
    models={"degree_days": linear_degree_days},
    start_date=start_date,
)
seq = pandas.date_range(start=start_date, end=end_date, freq="H")
TTmodel = DegreeDayModel(Tbase=0)
every_dd = thermal_time_filter(seq, weather, TTmodel, delay=20)
canopy_timing = IterWithDelays(*time_control(seq, every_dd, weather.data))

# Prepare saving of outputs
recorders = initiate_all_adel_septo_recorders(g, nsect)

# Reconstruction
for canopy_iter in canopy_timing:
    if canopy_iter:
        date = canopy_iter.value.index[-1]
        print(date)

        # Grow wheat
        g = adel.grow(g, canopy_iter.value)

        # Temporary : 3D display
        # scene = plot3d(g)
        # Viewer.display(scene)

        # Save variables
        for plant in recorders:
            for lf, recorder in recorders[plant].items():
                recorder.update_vids_with_labels(adel_ids=adel_ids(g))
                recorder.record_only_leaf_data(
                    g, date, degree_days=canopy_iter.value.degree_days[-1]
                )

scene = plot3d(g)
Viewer.display(scene)

for plant in recorders:
    for recorder in recorders[plant].values():
        recorder.create_dataframe_only_leaf_data()

leaves = ["F%d" % lf for lf in range(13)]
fig, axs = plt.subplots(3, 3)
for plant in recorders:
    for recorder in recorders[plant].values():
        data = recorder.data
        axs[0][0].plot(data.degree_days, data.leaf_area)
        axs[0][0].set_ylabel("leaf area (cm2)", fontsize=18)
        axs[0][1].plot(data.degree_days, data.leaf_green_area)
        axs[0][1].set_ylabel("leaf green area (cm2)", fontsize=18)
        axs[0][2].plot(data.degree_days, data.senesced_area)
        axs[0][2].set_ylabel("leaf senesced area (cm2)", fontsize=18)
        axs[1][0].plot(data.degree_days, data.leaf_length)
        axs[1][0].set_ylabel("leaf length (cm)", fontsize=18)
        axs[1][1].plot(data.degree_days, data.leaf_green_length)
        axs[1][1].set_ylabel("leaf green length (cm)", fontsize=18)
        axs[1][2].plot(data.degree_days, data.leaf_senesced_length)
        axs[1][2].set_ylabel("leaf senesced length (cm2)", fontsize=18)
        axs[2][0].plot(data.degree_days, data.leaf_min_height)
        axs[2][0].set_ylabel("leaf bottom height (cm)", fontsize=18)
        axs[2][1].plot(data.degree_days, data.leaf_mean_height)
        axs[2][1].set_ylabel("leaf mean height (cm)", fontsize=18)
        axs[2][2].plot(data.degree_days, data.leaf_max_height)
        axs[2][2].set_ylabel("leaf top height (cm)", fontsize=18)
        ylims = [
            min([ax.get_ylim()[0] for ax in axs[2]]),
            max([ax.get_ylim()[1] for ax in axs[2]]),
        ]
        for ax in axs[2]:
            ax.set_ylim(ylims)
