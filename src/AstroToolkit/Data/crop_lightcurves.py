def crop_lightcurve(struct, start=None, stop=None, timeformat=None):
    data = struct.data

    def filter_data(band, mask):
        return [val for i, val in enumerate(band) if i not in mask]

    combined_time = []
    for band in data:
        combined_time += band["mjd"]
    sync_time = min(combined_time)
    max_time = max(combined_time)

    if start[-1] == "%":
        start = float(start.rstrip("%"))
        start = sync_time + start / 100 * (max_time - sync_time)
    else:
        start = float(start)
        if timeformat == "reduced":
            start += sync_time
    if stop[-1] == "%":
        stop = float(stop.rstrip("%"))
        stop = sync_time + stop / 100 * (max_time - sync_time)
    else:
        stop = float(stop)
        if timeformat == "reduced":
            stop += sync_time

    """ Old
    if timeformat == "reduced" and kind == "time":
        start += sync_time
        stop += sync_time
    if kind == "percent":
        start = sync_time + start / 100 * (max_time - sync_time)
        stop = sync_time + stop / 100 * (max_time - sync_time)
    """

    cropped_data = []
    for i, band in enumerate(data):
        # skips any bands that did not return any data
        if not band["mag"]:
            return band

        # original time
        time = band["mjd"]

        bad_indices = [i for i, val in enumerate(time) if val < start or val > stop]

        # if no data within requested range in time, band is empty -> None
        data_remaining = True
        if len(bad_indices) == len(band["mag"]):
            data_remaining = False

        for key in band:
            if data_remaining:
                if key != "band":
                    band[key] = filter_data(band[key], bad_indices)
            else:
                if key != "band":
                    band[key] = None

        cropped_data.append(band)

    return cropped_data
