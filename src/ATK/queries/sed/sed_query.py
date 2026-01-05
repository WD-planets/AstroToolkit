import numpy as np
import pandas as pd

from ...structures.definitions import SED, Target
from ...Tools import query as general_query
from ...utilities.defaults import RETURNS
from .sed_core import SED_INFO, ab_mag_to_flux_mjy, get_ab_mag_offset


def get_survey_phot(survey: str, survey_data: pd.DataFrame) -> pd.DataFrame:
    """
    Return SED photometry for ALL VizieR detections of a survey, converts Vega -> AB magnitude (if needed) -> flux
    """

    info = SED_INFO[survey]

    # identify bands
    mag_cols = info["mag_names"]
    err_cols = info["err_names"]
    wavelengths = info["lambda_ref"]

    sed_rows = []
    for i, (mag_col, err_col, wl) in enumerate(zip(mag_cols, err_cols, wavelengths)):
        if mag_col not in survey_data.columns:
            continue

        if "_r" in survey_data:
            band_df = survey_data[["_r", mag_col, err_col]].copy()
        else:
            # Gaia source queries don't have _r since not using a cone search
            band_df = survey_data[[mag_col, err_col]].copy()
            band_df["_r"] = np.nan

        band_df = band_df.rename(columns={mag_col: "mag", err_col: "mag_err"})

        band_df["band"] = mag_col
        band_df["wavelength"] = wl
        band_df["survey"] = survey

        # Drop missing magnitudes
        band_df = band_df[np.isfinite(band_df["mag"])]

        # Vega -> AB if needed
        if "zp_vega" in info:
            ab_offset = get_ab_mag_offset(info["zp_vega"][i])
            band_df["mag_ab"] = band_df["mag"] + ab_offset
        else:
            band_df["mag_ab"] = band_df["mag"]

        # AB mag -> flux (mJy)
        band_df["flux_mjy"] = ab_mag_to_flux_mjy(band_df["mag_ab"])

        # flux err
        band_df["flux_err_mjy"] = band_df["flux_mjy"] * (np.log(10.0) / 2.5) * band_df["mag_err"]

        sed_rows.append(band_df)

    if not sed_rows:
        return pd.DataFrame()

    sed = pd.concat(sed_rows, ignore_index=True)

    return sed


def query(target: Target, radius: float, **kwargs):
    """
    Constructs an SED by combining photometry from Vizier catalogues
    """

    sed_tables = []

    # perform queries
    for survey in SED_INFO:
        data = general_query(kind="vizier", survey=survey, target=target, radius=radius)

        if data.exception or not data.data:
            return data

        # get SED dataframe for each survey
        phot = get_survey_phot(survey, data.data[0])
        if not phot.empty:
            sed_tables.append(phot)

    if not sed_tables:
        return RETURNS.NULL

    # combine surveys
    df = pd.concat(sed_tables, ignore_index=True)

    sed = SED(
        survey=df["survey"].to_numpy(),
        band=df["band"].to_numpy(),
        wavelength=df["wavelength"].to_numpy(),
        flux=df["flux_mjy"].to_numpy(),
        flux_err=df["flux_err_mjy"].to_numpy(),
        separation=df["_r"].to_numpy(),
    )

    return sed
