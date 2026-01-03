import pandas as pd

from ...structures.definitions import Lightcurve

# required light curve columns (in order)
REQUIRED_COLS = ["band", "mjd", "mag", "mag_err", "ra", "dec"]


def get_lightcurves(survey: str, data: pd.DataFrame) -> Lightcurve:
    # check for any missing columns
    missing_cols = [col for col in REQUIRED_COLS if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Lightcurve DataFrame missing required columns {', '.join(missing_cols)}.")

    # keep only required columns + sort
    data = data[REQUIRED_COLS]

    lcs = []
    for band in data["band"].unique():
        band_data = data[data.band == band]

        lcs.append(Lightcurve.from_dataframe(band_data, survey=survey, band=band))

    return lcs
