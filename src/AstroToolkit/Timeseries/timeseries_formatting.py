import pandas as pd


def format_data(data):
    data_exists = False
    for band in data:
        if band["mag"] is not None:
            data_exists = True

    if not data_exists:
        return None

    dataframes = [
        pd.DataFrame.from_dict(band) for band in data if band["mag"] is not None
    ]
    combined_data = pd.concat(dataframes)
    if "hjd_ori" in list(combined_data.keys()):
        combined_data.sort_values("hjd_ori", inplace=True)
        x_data = combined_data["hjd_ori"] - 2400000.5
    elif "mjd_ori" in list(combined_data.keys()):
        x_data = combined_data["mjd"]
    else:
        raise Exception("Failed to read input data structure.")

    combined_data.reset_index(drop=True, inplace=True)

    error_data = combined_data["mag_err"]

    filters = list(set(combined_data["band"]))
    first_filter = filters[0]
    med_first_filter = combined_data.query("band==@first_filter")["mag"].median()

    for filter in filters:
        med_filtered = combined_data.query("band==@filter")["mag"].median()
        combined_data.loc[combined_data["band"] == filter, "mag"] += (
            med_first_filter - med_filtered
        )

    y_data = combined_data["mag"]

    return list(x_data), list(y_data), list(error_data)
