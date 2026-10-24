import numpy as np
import pandas
from pandas.core.datetools import Hour
import math
import datetime


def thermal_time(
    requested_dates,
    emergence_date,
    input_data,
    data_type,
    thermal_time_increment_method,
    latitude,
):
    """Calculates the thermal time increment.

    :Parameters:
        - `requested_dates` : The dates for which the thermal time increment has to be calculated. These dates must belong to input_data, and be >= to emergence_date.
        - `emergence_date` : The emergence date. The interval [emergence date, emergence date + 1 year[ must be included in input_data. emergence_datetime must be <= than requested_dates
        - `input_data` : The input data. Expect a pandas.DataFrame object. The dataframe index must contain dates.
        - `data_type` : The type of the input data. Can be either 'daily' or 'hourly'.
        - `thermal_time_increment_method` : The method used for thermal time increment calculation. Can be either 'linear_TT' or 'compensated_TT'.
        - `latitude` : The latitude where the data have been measured. Used to convert 'daily' data to 'hourly' data.

    :Types:
        - `requested_dates` : pandas.DataFrame
        - `emergence_date` : dict
        - `input_data` : pandas.DataFrame
        - `data_type` : str enumerate
        - `thermal_time_increment_method` : str
        - `latitude` : float

    :return: The thermal time increment calculated.
    :rtype: pandas.DataFrame

    """
    print("emergence_date", emergence_date)
    emergence_datetime = pandas.datetime(
        emergence_date["year"],
        emergence_date["month"],
        emergence_date["day"],
        emergence_date["hour"],
        emergence_date["minute"],
        emergence_date["second"],
    )
    # check that input_data contains 'emergence_datetime' and 'emergence_datetime + 1 year'
    if emergence_datetime not in input_data.index:
        raise Exception(
            " ".join(
                [
                    "emergence_datetime",
                    emergence_datetime,
                    "does not belong to input_data",
                ]
            )
        )
    start = emergence_datetime
    end = pandas.datetime(
        emergence_datetime.year + 1,
        emergence_datetime.month,
        emergence_datetime.day,
        emergence_datetime.hour,
        emergence_datetime.minute,
        emergence_datetime.second,
    )
    if end not in input_data.index:
        raise Exception(
            " ".join(
                ["emergence_datetime + 1 year", end, "does not belong to input_data"]
            )
        )

    def crit_year(index_i):
        # Permits to extract data from start to end
        if start <= index_i <= end:
            return True
        return False

    # extract data from start to end
    extracted_data = input_data.select(crit_year)

    if data_type == "daily":
        # transform to hourly data
        hourly_data = parton_logan(extracted_data, latitude)
    else:  # hourly data
        # calculate the mean for each (Tmin,Tmax) tuple
        Tair_mean = (extracted_data["Tmin"] + extracted_data["Tmax"]) / 2.0
        hourly_data = pandas.DataFrame(
            Tair_mean, index=extracted_data.index, columns=["Tair"]
        )

    # check that input_data contains requested_dates
    # and that emergence_datetime is <= to requested_dates
    for requested_date in requested_dates.index:
        print("requested_date", type(requested_date), requested_date)
        if requested_date not in input_data.index:
            raise Exception(
                " ".join(
                    ["requested_date", requested_date, "does not belong to input_data"]
                )
            )
        if requested_date < emergence_datetime:
            raise Exception(
                " ".join(
                    [
                        "requested_date",
                        requested_date,
                        "is < than emergence_datetime",
                        emergence_datetime,
                    ]
                )
            )

    # calculate thermal time increment for all dates with the given thermal_time_increment_method
    available_methods = {"linear_TT": linear_TT, "compensated_TT": compensated_TT}
    thermal_time_increment_array = available_methods[thermal_time_increment_method](
        hourly_data["Tair"].values
    )
    thermal_time_increment_dataframe = pandas.DataFrame(
        thermal_time_increment_array,
        index=hourly_data.index,
        columns=[thermal_time_increment_method],
    )

    def crit_dates(index_i):
        # Permits to select requested_dates
        if index_i in requested_dates.index:
            return True
        return False

    # select data for requested_dates and return them
    return thermal_time_increment_dataframe.select(crit_dates)


def parton_logan(dailyMinMaxTemp, latitude=55, param="air150cm"):
    """Estimate hourly temperature from daily temperature. Based on Parton & Logan formula (AFM 1981)

    :Parameters:
        - `dailyMinMaxTemp` : the min and max temperatures for each consecutive julian day number.
        - `latitude` : latitude in degrees
        - `param` : a dictionary of 3 coefficients 'a', 'b' and 'c' (e.g.: {'a': 1.86, 'b': 2.2, 'c': - 0.17})
        OR a string indicating the type of temperature recorded ("air150cm", "air10cm", "soilsurface", or "soil-10cm").
        a,b,c are for :
        a = time lag coeff for max Temperature after noon (h)
        b = coefficient for temperature decrease at night
        c = timeLag for the minimum temperature after sunrise (h)

    :Types:
        - `dailyMinMaxTemp` : pandas.DataFrame
        - `latitude` : float
        - `param` : dictionary

    :return: the air temperature at hourly time steps.
    :rtype: pandas.DataFrame

    """

    paramref = {
        "air150cm": {"a": 1.86, "b": 2.2, "c": -0.17},
        "air10cm": {"a": 1.52, "b": 2, "c": -0.18},
        "soilsurface": {"a": 0.5, "b": 1.81, "c": 0.49},
        "soil-10cm": {"a": 0.45, "b": 2.28, "c": 1.83},
    }

    if type(param) is str:
        p = paramref[param]
    else:
        p = param

    a = p["a"]
    b = p["b"]
    c = p["c"]

    phi = latitude / 180.0 * math.pi

    Tmin = dailyMinMaxTemp["Tmin"]
    Tmax = dailyMinMaxTemp["Tmax"]

    daily_series = pandas.Series(list(range(Tmin.size)), index=Tmin.index)

    def calc_daylength(i):
        dayNumber = (
            dailyMinMaxTemp.index[i].toordinal()
            + 1
            - datetime.datetime(dailyMinMaxTemp.index[i].year, 1, 1).toordinal()
        )
        delta_i = 0.4014 * math.sin(2 * math.pi * (dayNumber - 77) / 365.0)
        t2_i = (-math.tan(phi)) * math.tan(delta_i)
        t1_i = math.sqrt(1 - t2_i**2)
        # angular amplitude
        ahr_i = math.atan2(t1_i, t2_i)
        return ahr_i / math.pi * 24

    daylength = daily_series.map(calc_daylength)

    def calc_sunrise(i):
        return 12 - daylength[i] / 2.0

    sunrise = daily_series.map(calc_sunrise)

    def calc_sunset(i):
        return 12 + daylength[i] / 2.0

    sunset = daily_series.map(calc_sunset)

    # minimal temperature hour
    def calc_hmin(i):
        return sunrise[i] + c

    hmin = daily_series.map(calc_hmin)

    # sunset temperature
    def calc_Tsunset(i):
        return Tmin[i] + (Tmax[i] - Tmin[i]) * math.sin(
            math.pi * (sunset[i] - hmin[i]) / (daylength[i] + 2 * a)
        )

    Tsunset = daily_series.map(calc_Tsunset)

    # sunset temperature at the day before
    Tsunsetb = pandas.Series(
        np.concatenate(([Tsunset.values[0]], Tsunset.values[:-1])), index=Tsunset.index
    )
    # minimal temperature at the day after
    Tmina = pandas.Series(
        np.concatenate((Tmin.values[1:], [Tmin.values[-1]])), index=Tsunset.index
    )

    hourly_idx = pandas.DateRange(
        Tmin.index[0], Tmin.index[-1] + 23 * Hour(), offset=Hour()
    )

    hourly_series = pandas.Series(list(range(hourly_idx.size)), index=hourly_idx)

    def calc_hourlyTemp(abs_hour):
        abs_day = abs_hour / 24
        rel_hour = abs_hour % 24 + 1
        if rel_hour < sunrise[abs_day]:
            return Tmin[abs_day] + (Tsunsetb[abs_day] - Tmin[abs_day]) * math.exp(
                -b * (24 - sunset[abs_day] + rel_hour) / (24 - daylength[abs_day])
            )
        elif rel_hour > sunset[abs_day]:
            return Tmina[abs_day] + (Tsunset[abs_day] - Tmina[abs_day]) * math.exp(
                -b * (rel_hour - sunset[abs_day]) / (24 - daylength[abs_day])
            )
        else:
            return Tmin[abs_day] + (Tmax[abs_day] - Tmin[abs_day]) * math.sin(
                math.pi * (rel_hour - hmin[abs_day]) / (daylength[abs_day] + 2 * a)
            )

    temp_series = hourly_series.map(calc_hourlyTemp)
    return pandas.DataFrame(temp_series, index=temp_series.index, columns=["Tair"])


def linear_TT(Th, Tb=0.0):
    """???

    :Parameters:
        - `Th` : ???
        - `Tb` : ???

    :Types:
        - `Th` : np.array
        - `Tb` : float

    :return: Cumulated hourly thermal time increment.
    :rtype: np.array

    """
    dsT_24 = np.maximum(np.zeros_like(Th), Th - Tb) / 24.0
    return dsT_24.cumsum()


def loiT(TK, k, EaR, DSR, DHR):
    """???

    :Parameters:
        - `TK` : temperature in Kelvin
        - `k` : ???
        - `EaR` : ???
        - `DSR` : ???
        - `DHR` : ???

    :Types:
        - `TK` : float
        - `k` : float
        - `EaR` : float
        - `DSR` : float
        - `DHR` : float

    :return: ???
    :rtype: float

    """
    return k * TK * np.power(math.e, -EaR / TK) / (1 + np.power(math.e, DSR - DHR / TK))


def compensated_TT(Tair, TCref=12.0, k=3.8088e10, EaR=8899.6, DSR=68.0, DHR=20736.0):
    """Hourly compensated thermal time increment.

    :Parameters:
        - `Tair` : temperature in Kelvin
        - `TCref` : ???
        - `k` : ???
        - `EaR` : ???
        - `DSR` : ???
        - `DSR` : ???

    :Types:
        - `Tair` : np.array
        - `TCref` : float
        - `k` : float
        - `EaR` : float
        - `DSR` : float
        - `DSR` : float

    :return: Cumulated hourly thermal time compensated increment.
    :rtype: np.array

    """
    TKref = 273.0 + TCref
    TKref_mod = loiT(TKref, k, EaR, DSR, DHR)
    TK_array = Tair + 273.0
    TK_array_mod = loiT(TK_array, k, EaR, DSR, DHR)
    dsTc = TK_array_mod * TCref / TKref_mod / 24.0
    return dsTc.cumsum()
