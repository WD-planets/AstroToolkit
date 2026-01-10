import astropy.units as u
import pandas as pd
from astropy.time import Time
from pyasassn.client import SkyPatrolClient

from ...structures.definitions import Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from .lightcurve_core import get_lightcurves


def query(target: Target, **kwargs: dict):
    client = SkyPatrolClient(verbose=False)

    # radius in deg
    radius = kwargs["radius"].to(u.deg).value

    try:
        data = client.cone_search(
            ra_deg=target.coords.ra.value, dec_deg=target.coords.dec.value, radius=radius, catalog="master_list", download=True
        )
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    ids = data.catalog_info["asas_sn_id"]

    all_lcs = []
    for id in ids:
        lc = data[id].data
        df = lc.rename(columns={"phot_filter": "band", "asas_sn_id": "id"})
        df["mjd"] = Time(lc["jd"], format="jd").mjd
        df["ra"] = data.catalog_info["ra_deg"][0]
        df["dec"] = data.catalog_info["dec_deg"][0]
        df["band"] = df["band"].str.lower()
        df["id"] = id

        # bad detections have mag_err = 100
        df = df[df["mag_err"] < 99]

        all_lcs.append(df)
    combined_lcs = pd.concat(all_lcs)

    lcs = get_lightcurves("asassn", target, kwargs["radius"], combined_lcs, kwargs.get("split", False))

    return lcs
