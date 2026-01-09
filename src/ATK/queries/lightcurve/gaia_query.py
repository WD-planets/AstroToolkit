import numpy as np
import pandas as pd
from astropy.time import Time

from ...structures.definitions import Target
from ...Tools import query as general_query
from ...utilities.defaults import RETURNS
from ..vizier.vizier_query import gaia_query_by_source
from .lightcurve_core import get_lightcurves


def query(target: Target, **kwargs: dict):
    """
    Performs a Gaia light curve query
    """

    if target.identifier and target.survey == "gaia":
        lc_data = gaia_query_by_source(target.identifier, kind="lightcurve")
        if lc_data is RETURNS.NULL or lc_data is RETURNS.EXCEPTION:
            return lc_data
    else:
        lc_data = general_query(kind="vizier", target=target, catalogue="I/355/epphot", radius=kwargs["radius"])
        # data returned
        if lc_data.data:
            lc_data = lc_data.data[0]
        # no data returned, exception encountered
        elif lc_data.exception:
            return RETURNS.EXCEPTION
        # no data returned, no exception encountered
        else:
            return RETURNS.NULL

    all_bands = []

    for band in ["G", "BP", "RP"]:
        band_data = pd.DataFrame(
            {
                "band": band.lower(),
                "ra": lc_data["RA_ICRS"],
                "dec": lc_data["DE_ICRS"],
                "mag": lc_data[f"{band}mag"],
                "time": lc_data[f"Time{band}"],
                "id": lc_data["Source"],
            }
        )

        # calculate magnitude errors from fluxes and errors
        flux = lc_data[f"F{band}"]
        flux_err = lc_data[f"e_F{band}"]
        mag_err = np.full(len(flux), np.nan)
        mask = (flux > 0) & (flux_err > 0)
        mag_err[mask] = (2.5 / np.log(10)) * (flux_err[mask] / flux[mask])
        band_data["mag_err"] = mag_err

        # calculate MJD from per-band time (need to remove nan times first)
        band_data = band_data.dropna(subset=["time"])
        band_data["mjd"] = Time(band_data["time"] + 2455197.5, format="jd").mjd

        all_bands.append(band_data)

    # combine bands into single DataFrame
    combined_data = pd.concat(all_bands)

    lcs = get_lightcurves(target, "gaia", combined_data, kwargs.get("split", False))

    return lcs
