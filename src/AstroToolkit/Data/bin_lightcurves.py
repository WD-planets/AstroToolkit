def binning(data, bins=None, bin_size=None):
    import numpy as np

    ra = data["ra"]
    dec = data["dec"]

    mag = data["mag"]
    mag_err = data["mag_err"]
    time = data["mjd"]

    sync_time = min(time)

    bin_edges_arr = []

    # if a number of equally spaced bins in time are requested
    if bins:
        bin_edge_lower = 0
        bin_size = (max(time) - min(time)) / bins
        bin_edge_higher = bin_size
        bins = np.linspace(0, bins, bins + 1)
        bin_edges = (bins * bin_size) + min(time)

        for i, val in enumerate(bin_edges):
            if i == len(bin_edges) - 1:
                break
            else:
                bin_edges_arr.append([bin_edges[i], bin_edges[i + 1]])

    # if a specific bin size is requested (e.g. '10h')
    elif bin_size:
        bin_edge_lower = min(time)

        bin_unit = bin_size[-1]
        bin_size = float(bin_size[:-1])

        # convert bin size to days based on unit provided
        if bin_unit == "d":
            bin_size = bin_size
        elif bin_unit == "h":
            bin_size = bin_size / 24
        elif bin_unit == "m":
            bin_size = bin_size / (24 * 60)

        bin_edge_higher = bin_edge_lower + bin_size

        # sets lower and higher bin edges using bin size in days, final one may be smaller than requested size if it exceeds upper time limit of data
        while True:
            reached_end = False
            if bin_edge_higher + bin_size >= max(time):
                reached_end = True

            bin_edges_arr.append([bin_edge_lower, bin_edge_higher])
            bin_edge_lower = bin_edge_higher
            bin_edge_higher += bin_size

            if reached_end:
                break

    # combine binned data
    final_mags, final_times, final_errors, final_ra, final_dec = [], [], [], [], []
    for bin_edges in bin_edges_arr:
        mask = [i for i, val in enumerate(time) if bin_edges[0] <= val < bin_edges[1]]

        if len(mask) == 0:
            continue

        binned_mags = [val for i, val in enumerate(mag) if i in mask]
        binned_ra = [val for i, val in enumerate(ra) if i in mask]
        binned_dec = [val for i, val in enumerate(dec) if i in mask]

        # calculates weighted mean based on errors (i.e. those with larger errors have less impact)
        weights = [1 / pow(err, 2) for i, err in enumerate(mag_err) if i in mask]
        weighted_mean, norm = np.average(binned_mags, weights=weights, returned=True)

        mean_ra = sum(binned_ra) / len(binned_ra)
        mean_dec = sum(binned_dec) / len(binned_dec)

        final_mags.append(weighted_mean)
        final_ra.append(mean_ra)
        final_dec.append(mean_dec)

        # sets time values to the middle of each bin
        final_times.append((bin_edges[1] + bin_edges[0]) / 2)
        # error on the weighted mean
        final_errors.append(1 / np.sqrt(norm))

    # set up final data dict
    data["mag"] = final_mags
    data["mjd"] = final_times
    data["mag_err"] = final_errors
    data["ra"] = final_ra
    data["dec"] = final_dec

    return data
