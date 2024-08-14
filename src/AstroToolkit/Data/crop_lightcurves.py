def crop_lightcurve(
    struct,
    start=None,
    stop=None,
    start_percent=None,
    stop_percent=None,
    timeformat="reduced",
):
    input_error = "Light curve cropping requires either: start, stop or start_percent, stop_percent"

    data = struct.data

    def filter_data(data, mask):
        return [val for i, val in enumerate(data) if i not in mask]

    def crop(data):
        cropped_data = []
        for band in data:
            if band["mag"] is None:
                cropped_data.append(band)
                continue

            if timeformat == "reduced":
                if "mjd" in band:
                    time = band["mjd"]
                elif "hjd" in band:
                    time = band["hjd"]
                else:
                    raise Exception("Could not find time data in given band.")
            if timeformat == "original":
                if "mjd_ori" in band:
                    time = band["mjd_ori"]
                elif "hjd_ori" in band:
                    time = band["hjd_ori"]
                else:
                    raise Exception("Could not find time data in given band.")

            bad_indices = [i for i, val in enumerate(time) if val < start or val > stop]

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

    def crop_percent(data):
        cropped_data = []
        for band in data:
            if band["mag"] is None:
                cropped_data.append(band)
                continue

            if timeformat == "reduced":
                if "mjd" in band:
                    time = band["mjd"]
                elif "hjd" in band:
                    time = band["hjd"]
                else:
                    raise Exception("Could not find time data in given band.")
            if timeformat == "original":
                if "mjd_ori" in band:
                    time = band["mjd_ori"]
                elif "hjd_ori" in band:
                    time = band["hjd_ori"]
                else:
                    raise Exception("Could not find time data in given band.")

            bad_indices = [
                i
                for i, val in enumerate(time)
                if val < (start_percent / 100 * max(time))
                or val > (stop_percent / 100 * max(time))
            ]

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

    if (start is not None and stop is not None) and (
        start_percent is None and stop_percent is None
    ):
        return crop(data)
    elif (start_percent is not None and stop_percent is not None) and (
        start is None and stop is None
    ):
        return crop_percent(data)
    else:
        raise ValueError(input_error)
