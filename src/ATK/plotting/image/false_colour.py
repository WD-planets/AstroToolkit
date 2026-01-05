import matplotlib as mpl
import numpy as np

N_COLOURS = 256

LAMBDA_MIN = 380
LAMBDA_MAX = 750
GAMMA = 0.8
OFFSET = 0.3
MULTIPLIER = 0.7

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


def wavelength_to_rgb(wavelength) -> tuple:
    """
    Performs an approximate conversion from wavelength in nm to sRGB
    """

    # clip extreme wavelengths
    wavelength = np.clip(wavelength, LAMBDA_MIN, LAMBDA_MAX)

    # (very roughly) convert wavelengths to RGB values
    if 380 <= wavelength <= 440:
        R = -(wavelength - 440.0) / (440.0 - 380.0)
        G = 0.0
        B = 1.0
    elif 440 < wavelength <= 490:
        R = 0.0
        G = (wavelength - 440.0) / (490.0 - 440.0)
        B = 1.0
    elif 490 < wavelength <= 510:
        R = 0.0
        G = 1.0
        B = -(wavelength - 510.0) / (510.0 - 490.0)
    elif 510 < wavelength <= 580:
        R = (wavelength - 510.0) / (580.0 - 510.0)
        G = 1.0
        B = 0.0
    elif 580 < wavelength <= 645:
        R = 1.0
        G = -(wavelength - 645.0) / (645.0 - 580.0)
        B = 0.0
    else:
        R = 1.0
        G = 0.0
        B = 0.0

    # intensity correction near vision limits
    if wavelength < 420:
        factor = OFFSET + MULTIPLIER * (wavelength - 380) / (420 - 380)
    elif wavelength > 645:
        factor = OFFSET + MULTIPLIER * (750 - wavelength) / (750 - 645)
    else:
        factor = 1.0

    # set gamma
    R = (R * factor) ** GAMMA
    G = (G * factor) ** GAMMA
    B = (B * factor) ** GAMMA

    return (R, G, B)


def get_false_cmap(survey: str, band: str) -> list:
    """
    Converts the effective wavelength of a given band of an imaging survey to a colour palette
    """

    wavelength = EFFECTIVE_WAVELENGTHS[survey][band]

    # convert wavelength to srgb
    colour = wavelength_to_rgb(wavelength)

    # convert to palette
    cmap = mpl.colors.LinearSegmentedColormap.from_list(f"{int(wavelength)}nm", [(0, "black"), (1, colour)])
    colours = [mpl.colors.rgb2hex(cmap(i / (N_COLOURS - 1))) for i in range(N_COLOURS)]

    return colours
