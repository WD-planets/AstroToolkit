def binning(data, bins=None, bin_size=None):
    import numpy as np

    ra = data["ra"]
    dec = data["dec"]

    mag = data["mag"]
    mag_err = data["mag_err"]
    if "mjd_ori" in data:
        time = data["mjd_ori"]
        time_unit = "mjd_ori"
    else:
        time = data["hjd_ori"]
        time_unit = "hjd_ori"

    bin_edges_arr = []

    if bins:
        # fixes a bug that made the script get stuck when this condition wasn't met
        if bins > len(mag):
            return data

        bin_edge_lower = 0
        bin_size = int(len(time) / bins)
        bin_edge_higher = bin_size

        while True:
            reached_end = False
            if bin_edge_higher + bin_size >= len(time):
                bin_edge_higher = len(time) - 1
                reached_end = True

            bin_edges_arr.append([bin_edge_lower, bin_edge_higher])
            bin_edge_lower = bin_edge_higher
            bin_edge_higher += bin_size

            if reached_end:
                break

        for i, val in enumerate(bin_edges_arr):
            bin_edges_arr[i][0] = time[val[0]]
            bin_edges_arr[i][1] = time[val[1]]

    elif bin_size:
        bin_edge_lower = min(time)

        try:
            bin_unit = bin_size[-1]
            bin_size = float(bin_size[:-1])
        except:
            raise Exception("No binsize unit given.")

        if bin_unit == "d":
            bin_size = bin_size
        elif bin_unit == "h":
            bin_size = bin_size / 24
        elif bin_unit == "m":
            bin_size = bin_size / (24 * 60)
        else:
            raise Exception("Invalid binsize unit.")

        bin_edge_higher = bin_edge_lower + bin_size

        while True:
            reached_end = False
            if bin_edge_higher + bin_size >= max(time):
                bin_edge_higher = max(time)
                reached_end = True

            bin_edges_arr.append([bin_edge_lower, bin_edge_higher])
            bin_edge_lower = bin_edge_higher
            bin_edge_higher += bin_size

            if reached_end:
                break

    final_mags, final_times, final_errors, final_ra, final_dec = [], [], [], [], []
    for bin_edges in bin_edges_arr:
        mask = [i for i, val in enumerate(time) if bin_edges[0] <= val < bin_edges[1]]

        if len(mask) == 0:
            continue

        binned_mags = [val for i, val in enumerate(mag) if i in mask]
        binned_ra = [val for i, val in enumerate(ra) if i in mask]
        binned_dec = [val for i, val in enumerate(dec) if i in mask]

        weights = [1 / pow(err, 2) for i, err in enumerate(mag_err) if i in mask]

        weighted_mean, norm = np.average(binned_mags, weights=weights, returned=True)

        mean_ra = sum(binned_ra) / len(binned_ra)
        mean_dec = sum(binned_dec) / len(binned_dec)

        final_mags.append(weighted_mean)
        final_ra.append(mean_ra)
        final_dec.append(mean_dec)
        final_times.append((bin_edges[1] + bin_edges[0]) / 2)
        final_errors.append(1 / np.sqrt(norm))

    base_time_unit = time_unit[:-4]

    data["mag"] = final_mags
    data[time_unit] = final_times
    data["mag_err"] = final_errors
    data[base_time_unit] = [t - min(final_times) for t in final_times]
    data["ra"] = final_ra
    data["dec"] = final_dec

    return data
