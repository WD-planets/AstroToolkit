import matplotlib as mpl

N_COLOURS = 256


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


def wavelength_to_cmap(wavelength_nm):
    colour = wavelength_to_rgb(wavelength_nm)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(f"{int(wavelength_nm)}nm", [(0, "black"), (1, colour)])
    colours = [mpl.colors.rgb2hex(cmap(i / (N_COLOURS - 1))) for i in range(N_COLOURS)]

    return colours
