import warnings
from io import BytesIO

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.io.fits import Header
from astropy.table import Table
from astropy.units import UnitsWarning
from astroquery.mast import Observations

from ...structures.Target import Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from ...utilities.misc import suppress_stdout
from ...utilities.requests import send_request
from .lightcurve_core import get_lightcurves

# ignore units warning when reading fits table
warnings.simplefilter("ignore", category=UnitsWarning)


def read_tess_data(data: pd.DataFrame, header: Header) -> pd.DataFrame:
    # create quality mask
    snr = data["PDCSAP_FLUX"] / data["PDCSAP_FLUX_ERR"]
    mask = np.isfinite(data["TIME"]) & np.isfinite(data["PDCSAP_FLUX"]) & np.isfinite(data["PDCSAP_FLUX_ERR"]) & (data["PDCSAP_FLUX"] > 0) & (snr > 3) & (data["QUALITY"] == 0)
    data = data[mask].copy()

    if not np.any(mask):
        return pd.DataFrame()

    # convert BTJD → MJD
    data["TIME"] += 2457000.0 - 2400000.5

    # calculate mag and mag_err
    data["mag"] = -2.5 * np.log10(data["PDCSAP_FLUX"]) + 20.44
    # 2.5/ln(10)
    data["mag_err"] = 1.085736 * (data["PDCSAP_FLUX_ERR"] / data["PDCSAP_FLUX"])

    # add band, ra, and dec columns
    data["band"] = "Tmag"
    data["ra"] = header.get("RA_OBJ")
    data["dec"] = header.get("DEC_OBJ")
    data["id"] = header.get("TICID")

    data = data.rename(columns={"TIME": "mjd"})

    return data


def query(target: Target, **kwargs: dict):
    """
    Performs a TESS light curve query
    """

    base_url = "https://mast.stsci.edu/api/v0.1/Download/file?uri="

    radius = kwargs["radius"].to(u.deg)

    # query region around target
    try:
        obs = Observations.query_region(target.coords, radius=radius)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION
    if not len(obs):
        return RETURNS.NULL

    # keep only TESS observations
    obs = obs[obs["obs_collection"] == "TESS"]
    if not len(obs):
        return RETURNS.NULL

    # get all products for found objects, need to suppress a logger-level warning here and some additional output
    with suppress_stdout():
        products_list = Observations.get_unique_product_list(obs)

    # get fits light curves for found objects
    lc_product_list = Observations.filter_products(products_list, productSubGroupDescription="LC", extension="fits")
    if not len(lc_product_list):
        return RETURNS.NULL

    rows = []
    for row in lc_product_list:
        # get mast URI
        mast_uri = row.get("dataURI")
        if not mast_uri:
            continue

        # construct URL and send request
        url = f"{base_url}{mast_uri}"
        response = send_request("tess", url)
        if response is RETURNS.EXCEPTION:
            continue

        # read response
        hdul = fits.open(BytesIO(response.content))
        lc_data = Table.read(hdul[1]).to_pandas()
        lc_header = hdul[0].header

        df = read_tess_data(lc_data, lc_header)

        rows.append(df)

    if not rows:
        return None

    df = pd.concat(rows, ignore_index=True)

    lcs = get_lightcurves("tess", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
