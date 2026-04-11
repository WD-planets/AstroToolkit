import warnings

warnings.filterwarnings("ignore", message=".*tpfmodel submodule is not available.*", category=UserWarning)

import logging
from contextlib import contextmanager

import astropy.units as u
import lightkurve as lk
import numpy as np
import pandas as pd
from astropy.io.fits import Header
from astropy.time import Time
from astropy.units import UnitsWarning
from lightkurve.search import SearchError

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from ...utilities.misc import suppress_stdout
from .lightcurve_core import get_lightcurves

# ignore units warning when reading fits table
warnings.simplefilter("ignore", category=UnitsWarning)


@contextmanager
def suppress_lightkurve():
    logger = logging.getLogger("lightkurve")
    default_level = logger.level

    # silence everything below CRITICAL
    logger.setLevel(logging.CRITICAL)

    # run code inside context manager
    try:
        yield

    # set back to default
    finally:
        logger.setLevel(default_level)


def read_tess_data(data: pd.DataFrame, header: Header, kwargs: dict) -> pd.DataFrame:
    # create quality mask
    flux, flux_err = data["flux"], data["flux_err"]

    snr = flux / flux_err

    # necessary
    mask = np.isfinite(data["time"]) & np.isfinite(flux) & np.isfinite(flux_err) & (flux > 0)

    # optional
    additional_mask = (snr > 3) & (data["quality"] == 0)
    final_mask = mask & additional_mask if kwargs.get("filter") else mask
    data = data[final_mask].copy()

    if not np.any(final_mask):
        return pd.DataFrame()

    # calculate mag and mag_err
    data["mag"] = -2.5 * np.log10(flux) + 20.44
    data["mag_err"] = 1.085736 * (flux_err / flux)

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

    try:
        with suppress_lightkurve():
            search = lk.search_lightcurve(target.coords.to_string("hmsdms"), mission="TESS", radius=radius)
    except Exception:
        return RETURNS.EXCEPTION

    if not len(search):
        return RETURNS.NULL

    preferred = ["TESS-SPOC", "QLP"]
    for p in preferred:
        subset = search[search.author == p]
        if len(subset) > 0:
            search = subset
            break

    if len(search) == 0:
        return RETURNS.NULL

    lcs = search.download_all()

    if lcs is None or len(lcs) == 0:
        return RETURNS.NULL

    rows = []

    for lc in lcs:
        # convert to table to pandas
        df = lc.to_table().to_pandas()

        df["time"] = lc.time.mjd
        df["flux"] = np.asarray(lc.flux, dtype=np.float64)
        df["flux_err"] = np.asarray(lc.flux_err, dtype=np.float64)
        header = lc.meta

        try:
            df = read_tess_data(df, header, kwargs)
        except Exception:
            continue

        if not df.empty:
            rows.append(df)

    if not rows:
        return RETURNS.NULL

    df = pd.concat(rows, ignore_index=True)

    lcs = get_lightcurves("tess", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
