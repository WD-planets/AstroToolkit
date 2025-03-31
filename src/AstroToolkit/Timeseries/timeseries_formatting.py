import pandas as pd

from ..Utility import getBrightnessType


def format_data(data):
    brightness_type = getBrightnessType(data)

    data_exists = False
    for band in data:
        if band[brightness_type]:
            data_exists = True

    if not data_exists:
        return None

    dataframes = [pd.DataFrame.from_dict(band) for band in data if band[brightness_type] is not None]
    combined_data = pd.concat(dataframes)
    x_data = combined_data["mjd"]

    combined_data.reset_index(drop=True, inplace=True)

    error_data = combined_data[f"{brightness_type}_err"]

    filters = list(set(combined_data["band"]))
    first_filter = filters[0]
    med_first_filter = combined_data.query("band==@first_filter")[brightness_type].median()

    for filter in filters:
        med_filtered = combined_data.query("band==@filter")[brightness_type].median()
        combined_data.loc[combined_data["band"] == filter, brightness_type] += med_first_filter - med_filtered

    y_data = combined_data[brightness_type]

    return list(x_data), list(y_data), list(error_data)
