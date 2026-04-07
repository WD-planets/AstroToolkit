import re

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.units import Quantity
from astroquery.simbad import Simbad

from ...utilities.defaults import CONNECTION_ERRORS, RETURNS


def get_ids(coords: list[SkyCoord], radius: Quantity):
    simbad = Simbad()
    simbad.ROW_LIMIT = -1

    if isinstance(coords, SkyCoord):
        targets_list = [coords]
    else:
        targets_list = list(coords)

    target_coords = SkyCoord(ra=[t.ra for t in targets_list], dec=[t.dec for t in targets_list], frame="icrs")

    try:
        response = simbad.query_region(target_coords, radius=radius)
    except Exception:
        return RETURNS.EXCEPTION

    if response is None or len(response) == 0:
        return RETURNS.NULL

    df = response.to_pandas()

    # clean IDs
    ids = np.asarray([re.sub(r"\s+", " ", str(x)) for x in df["main_id"]])

    simbad_coords = SkyCoord(ra=df["ra"].to_numpy() * u.deg, dec=df["dec"].to_numpy() * u.deg, frame="icrs")

    target_idx, simbad_idx, _, _ = search_around_sky(target_coords, simbad_coords, radius)

    matches = [[] for _ in range(len(targets_list))]

    for t_idx, s_idx in zip(target_idx, simbad_idx):
        matches[t_idx].append(ids[s_idx])

    # remove duplicates while preserving order
    matches = [list(dict.fromkeys(m)) for m in matches]

    # format output
    string_results = [", ".join(m) if m else "None" for m in matches]

    max_len = max(len(s) for s in string_results)
    results = np.asarray(string_results, dtype=f"<U{max_len}")

    return results
