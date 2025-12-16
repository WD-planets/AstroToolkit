import matplotlib as mpl

N_COLOURS = 256

# this will need updating, also haven't checked
EFFECTIVE_WAVELENGTHS = {
    "panstarrs": {"g": 481.0, "r": 615.5, "i": 750.3, "z": 866.8, "y": 961.4},
    "galex": {"fuv": 154.9, "nuv": 230.3},
    "sdss": {"u": 355.1, "g": 468.6, "r": 616.6, "i": 748.0, "z": 893.2},
    "2mass": {"j": 1069.1, "h": 1446.5, "ks": 2155.8},
    "wise": {"w1": 3400.0, "w2": 4600.0, "w3": 12000.0, "w4": 22000.0},
    "dss1": {"blue": 405.0, "red": 645.0},
    "dss2": {"blue": 480.0, "red": 670.0, "ir": 875.0},
    "skymapper": {"u": 350.0, "v": 450.0, "g": 480.0, "r": 625.0, "i": 775.0, "z": 870.0},
}


def wavelength_to_rgb(wavelength_nm):
    """
    Approximate conversion from wavelength to sRGB tuple.
    """

    if wavelength_nm < 380:
        wavelength_nm = 380
    if wavelength_nm > 750:
        wavelength_nm = 750
    gamma = 0.8

    if 380 <= wavelength_nm <= 440:
        R = -(wavelength_nm - 440.0) / (440.0 - 380.0)
        G = 0.0
        B = 1.0
    elif 440 < wavelength_nm <= 490:
        R = 0.0
        G = (wavelength_nm - 440.0) / (490.0 - 440.0)
        B = 1.0
    elif 490 < wavelength_nm <= 510:
        R = 0.0
        G = 1.0
        B = -(wavelength_nm - 510.0) / (510.0 - 490.0)
    elif 510 < wavelength_nm <= 580:
        R = (wavelength_nm - 510.0) / (580.0 - 510.0)
        G = 1.0
        B = 0.0
    elif 580 < wavelength_nm <= 645:
        R = 1.0
        G = -(wavelength_nm - 645.0) / (645.0 - 580.0)
        B = 0.0
    else:
        R = 1.0
        G = 0.0
        B = 0.0

    # intensity correction near vision limits
    if wavelength_nm < 420:
        factor = 0.3 + 0.7 * (wavelength_nm - 380) / (420 - 380)
    elif wavelength_nm > 645:
        factor = 0.3 + 0.7 * (750 - wavelength_nm) / (750 - 645)
    else:
        factor = 1.0

    R = (R * factor) ** gamma
    G = (G * factor) ** gamma
    B = (B * factor) ** gamma

    return (R, G, B)


def get_false_cmap(survey: str, band: str):
    wavelength_nm = EFFECTIVE_WAVELENGTHS[survey][band]

    colour = wavelength_to_rgb(wavelength_nm)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(f"{int(wavelength_nm)}nm", [(0, "black"), (1, colour)])
    colours = [mpl.colors.rgb2hex(cmap(i / (N_COLOURS - 1))) for i in range(N_COLOURS)]

    return colours
