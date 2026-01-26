from io import BytesIO

import astropy.units as u
import pandas as pd

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request
from .lightcurve_core import get_lightcurves


def query(target: Target, **kwargs: dict):
    """
    Performs a ZTF light curve query
    """

    # set up URL
    radius = kwargs["radius"].to(u.deg).value

    url = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE {target.coords.ra.value} {target.coords.dec.value} {radius}&BANDNAME=g,r,i&FORMAT=CSV"

    # get data
    response = send_request("ztf", url)
    if response is RETURNS.EXCEPTION:
        return response

    data = pd.read_csv(BytesIO(response.content))
    if not len(data):
        return RETURNS.NULL

    df = pd.DataFrame(
        {
            "mjd": data["mjd"],
            "mag": data["mag"],
            "mag_err": data["magerr"],
            "ra": data["ra"],
            "dec": data["dec"],
            "band": data["filtercode"].str[1:],
            "id": data["oid"],
        }
    )

    lcs = get_lightcurves("ztf", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
