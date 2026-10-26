"""
The user tutorial 05's example observations: SYNTHETIC, not measurements. They stand for what an experiment would give,
so that the tutorial can show how to give the model measured inputs, compare it with observations and fit a parameter.

- air_temperature.csv: a daily mean air temperature (°C), warming from 8 to 16 °C over 90 days, with day-to-day noise;
- observed_development.csv: every 10 days, the main stem's Haun stage, the plant's green leaf area (cm2) and its
  height (cm), made by the model at that temperature from the measured main stem shipped with Adel (data/
  Mariem_dimT_user_MIN.csv, plantgen's example input) with a phyllochron of 100 °Cd (the same example's dynT_user
  table: 0.010 leaves per °Cd), with noise (0.15 on the Haun stage, 5 % on the area and the height).

    python make_synthetic_data.py          (rewrites the two files in this folder)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

from openalea.adel.simulation import simulate, tables_from_measurements
from openalea.adel.tables import DATA

HERE = Path(__file__).parent
TRUE_PHYLLOCHRON = 100.
DAYS = 90


def main():
    noise = np.random.default_rng(0)
    days = np.arange(DAYS + 1)
    temperature = pd.DataFrame({"day": days, "temperature": np.round(
        8. + 8. * days / DAYS + noise.normal(0., 1.5, days.size), 1)})
    temperature.to_csv(HERE / "air_temperature.csv", index=False)

    dimensions = pd.read_csv(os.path.join(DATA, "Mariem_dimT_user_MIN.csv"))
    result = simulate(dict(tables_from_measurements(dimensions), phyllochron=TRUE_PHYLLOCHRON), days=DAYS,
                      temperature=temperature.set_index("day")["temperature"])
    observed = result.summary[(result.summary["day"] % 10 == 0) & (result.summary["day"] > 0)]
    observed = pd.DataFrame({
        "day": observed["day"].astype(int),
        "haun_stage": np.round(observed["haun_stage"] + noise.normal(0., 0.15, len(observed)), 2),
        "green_leaf_area": np.round(observed["green_leaf_area"] * 1e4 * noise.normal(1., 0.05, len(observed)), 1),
        "height": np.round(observed["height"] * 100 * noise.normal(1., 0.05, len(observed)), 1)})
    observed.to_csv(HERE / "observed_development.csv", index=False)


if __name__ == "__main__":
    main()
