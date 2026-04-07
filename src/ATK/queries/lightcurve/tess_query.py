import warnings

warnings.filterwarnings("ignore", message=".*tpfmodel submodule is not available.*", category=UserWarning)

import astropy.units as u
import lightkurve as lk
import numpy as np
import pandas as pd
from astropy.io.fits import Header
from astropy.time import Time
from astropy.units import UnitsWarning

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from .lightcurve_core import get_lightcurves

# ignore units warning when reading fits table
warnings.simplefilter("ignore", category=UnitsWarning)


def read_tess_data(data: pd.DataFrame, header: Header, kwargs: dict) -> pd.DataFrame:
    # create quality mask
    snr = data["pdcsap_flux"] / data["pdcsap_flux_err"]

    # necessary
    mask = np.isfinite(data["time"]) & np.isfinite(data["pdcsap_flux"]) & np.isfinite(data["pdcsap_flux_err"]) & (data["pdcsap_flux"] > 0)

    # optional
    additional_mask = (snr > 3) & (data["quality"] == 0)
    final_mask = mask & additional_mask if kwargs.get("filter") else mask
    data = data[final_mask].copy()

    if not np.any(final_mask):
        return pd.DataFrame()

    # convert BTJD → MJD
    data["time"] = Time(data["time"]).mjd.astype("float")

    # calculate mag and mag_err
    data["mag"] = -2.5 * np.log10(data["pdcsap_flux"]) + 20.44
    data["mag_err"] = 1.085736 * (data["pdcsap_flux_err"] / data["pdcsap_flux"])

    # add band, ra, and dec columns
    data["band"] = "Tmag"
    data["ra"] = header.get("RA_OBJ")
    data["dec"] = header.get("DEC_OBJ")
    data["id"] = header.get("TICID")

    data = data.rename(columns={"time": "mjd"})

    return data


def query(target: Target, **kwargs: dict):
    """
    Performs a TESS light curve query using lightkurve
    """

    radius = kwargs["radius"].to(u.deg)

    # ======================
    # SEARCH (lightkurve)
    # ======================

    try:
        search = lk.search_lightcurve(target.coords.to_string("hmsdms"), mission="TESS", radius=radius)
    except Exception:
        return RETURNS.EXCEPTION

    search = search[search.author == "TESS-SPOC"]

    if len(search) == 0:
        return RETURNS.NULL

    lcs = search.download_all()

    if lcs is None or len(lcs) == 0:
        return RETURNS.NULL

    rows = []

    for lc in lcs:
        # convert to table to pandas
        tab = lc.to_table()[["time", "flux", "flux_err", "pdcsap_flux", "pdcsap_flux_err", "quality"]]
        df = tab.to_pandas()

        header = lc.meta

        # ensure required columns exist
        if not all(col in df.columns for col in ["time", "pdcsap_flux", "pdcsap_flux_err", "quality"]):
            continue

        df = read_tess_data(df, header, kwargs)

        if not df.empty:
            rows.append(df)

    if not rows:
        return None

    df = pd.concat(rows, ignore_index=True)

    lcs = get_lightcurves("tess", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
