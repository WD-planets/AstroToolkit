import pandas as pd

from ...structures.definitions import Lightcurve, Target

# required light curve columns (in order)
REQUIRED_COLS = ["band", "mjd", "mag", "mag_err", "ra", "dec"]


def get_lightcurves(target: Target, survey: str, data: pd.DataFrame, split: bool) -> Lightcurve:
    # check for any missing columns
    missing_cols = [col for col in REQUIRED_COLS if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Lightcurve DataFrame missing required columns {', '.join(missing_cols)}.")

    if "id" not in data:
        split = False

    # split dataframe into light curves per-band
    lcs = []
    for band in data["band"].unique():
        band_data = data[data.band == band]
        if band_data.empty:
            continue

        # split light curves object ID if requested and id column exists
        if split:
            for id in data["id"].unique():
                obj_data = band_data[band_data.id == id]
                if obj_data.empty:
                    continue

                obj_data = obj_data[REQUIRED_COLS]
                lc = Lightcurve.from_dataframe(obj_data, survey=survey, band=band)
                lc.obj_id = str(id)
                lcs.append(lc)
        else:
            band_data = band_data[REQUIRED_COLS]
            lcs.append(Lightcurve.from_dataframe(band_data, survey=survey, band=band))

    # sort by object ID
    if split:
        lcs.sort(key=lambda lc: lc.obj_id)

    return lcs
