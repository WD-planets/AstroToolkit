import numpy as np
import pandas as pd

from AstroToolkit.Tools import query


def get_estimate(bprp, abs_g):
    suffix = "V"

    if 0.4 < bprp <= 0.8:
        prefix = "F"
        start = 0.4
        interval = 0.4
        num = 10
    elif 0.8 < bprp <= 1.0:
        prefix = "G"
        start = 0.8
        interval = 0.2
        num = 10
    elif 1.0 < bprp <= 1.85:
        prefix = "K"
        start = 1.0
        interval = 0.85
        num = 10
    elif 1.85 < bprp <= 2.95:
        prefix = "M"
        start = 1.85
        interval = 1.1
        num = 5
    else:
        raise Exception("Unexpected Error.")

    if 1.25 < bprp < 1.55 and 3.7 < abs_g < 4.6:
        suffix = "IV"
    elif 1.05 < bprp < 1.25 and 2.5 < abs_g < 3.6:
        suffix = "III"

    remainder = str(int((bprp - start) / (interval / num)))
    sType = f"{prefix}{remainder}{suffix}"
    return sType


def get_spectral_types(sources):
    for source in sources:
        gaia_data = query(kind="data", source=source, survey="gaia").data

        abs_g = (
            gaia_data["phot_g_mean_mag"][0]
            + 5 * np.log10(gaia_data["parallax"][0] / 1000)
            + 5
        )

    spectral_types = []
    spectral_types.append(
        get_estimate(
            gaia_data["phot_bp_mean_mag"][0] - gaia_data["phot_rp_mean_mag"][0], abs_g
        )
    )

    return spectral_types
