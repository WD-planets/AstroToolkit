"""
The function below (rolling_window_sigma_clip) is adapted from the source code for the fundamentals package: https://github.com/thespacedoctor/fundamentals/tree/master, however part of the package that ATK won't use
doesn't seem to work on windows (as it needs the readline package which isn't available for windows). This lets the package bypass the installation and still work on windows/linux
simultaneously.
"""


def rolling_window_sigma_clip(array, clippingSigma, windowSize):
    from astropy.stats import mad_std, sigma_clip

    midWindow = int((windowSize + 1) / 2)
    # ACCOMODATE SMALL LIST SIZES
    if len(array) < 5:
        return len(array) * [False]
    elif len(array) < windowSize:
        masked = sigma_clip(
            array,
            sigma_lower=clippingSigma,
            sigma_upper=clippingSigma,
            maxiters=7,
            cenfunc="median",
            stdfunc=mad_std,
        )
        return list(masked.mask)

    startOfWindow = 0
    endOfWindow = windowSize
    maskedArray = []
    dataIndex = 0
    while len(array) >= endOfWindow:
        arrayWindow = array[startOfWindow:endOfWindow]
        startOfWindow += 1
        endOfWindow += 1
        masked = sigma_clip(
            arrayWindow,
            sigma_lower=clippingSigma,
            sigma_upper=clippingSigma,
            maxiters=7,
            cenfunc="median",
            stdfunc=mad_std,
        )

        if dataIndex == 0:
            # 0,1,2...midWindow-1
            maskedArray += list(masked.mask)[0:midWindow]
            dataIndex += midWindow
        elif len(array) == endOfWindow - 1:
            # -midWindow...-2,-1
            maskedArray += list(masked.mask)[-midWindow:]
            dataIndex += midWindow
        else:
            maskedArray += [list(masked.mask)[midWindow - 1]]
            dataIndex += 1

    return maskedArray


def sigma_clip(data, sigma, window_size=11):
    from operator import itemgetter

    # sigma clip data with a rolling window
    if "hjd" in list(data.keys()):
        time_unit = "hjd"
    elif "mjd" in list(data.keys()):
        time_unit = "mjd"

    time_arr = data[time_unit]
    time_ori_arr = data[time_unit + "_ori"]
    mag_arr = data["mag"]
    mag_err_arr = data["mag_err"]
    ra_arr = data["ra"]
    dec_arr = data["dec"]

    # create list of dicts containing single data points
    sigma_clip_data = []
    for i in range(0, len(data["mag"])):
        sigma_clip_data.append(
            {
                time_unit: time_arr[i],
                f"{time_unit}_ori": time_ori_arr[i],
                "mag": mag_arr[i],
                "mag_err": mag_err_arr[i],
                "ra": ra_arr[i],
                "dec": dec_arr[i],
            }
        )

    # sort by mjd_ori
    sigma_clip_data = sorted(
        sigma_clip_data, key=itemgetter(f"{time_unit}_ori"), reverse=False
    )

    sigma_clip_mag = []
    sigma_clip_mag[:] = [row["mag"] for row in sigma_clip_data]

    fullMask = rolling_window_sigma_clip(
        array=sigma_clip_mag, clippingSigma=sigma, windowSize=window_size
    )

    try:
        sigma_clip_data = [e for e, m in zip(sigma_clip_data, fullMask) if not m]
    except:
        sigma_clip_data = []

    data['mag'] = [x["mag"] for x in sigma_clip_data]
    data['mag_err'] = [x["mag_err"] for x in sigma_clip_data]
    data[time_unit] = [x[time_unit] for x in sigma_clip_data]
    data[f'{time_unit}_ori'] = [x[f"{time_unit}_ori"] for x in sigma_clip_data]
    data['ra'] = [x["ra"] for x in sigma_clip_data]
    data['dec'] = [x["dec"] for x in sigma_clip_data]

    return data
