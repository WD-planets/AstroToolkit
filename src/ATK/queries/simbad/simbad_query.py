import re

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.units import Quantity
from astroquery.simbad import Simbad

from ...utilities.defaults import CONNECTION_ERRORS, RETURNS


def get_ids(targets: SkyCoord, radius: Quantity):
    """
    Fetch SIMBAD IDs for each detection in `targets`.

    Returns
    -------
    np.ndarray of str or None
        One element per target. Each element is a human-readable, list-like
        string of all SIMBAD IDs within `radius`, e.g.
        '["HD 12345", "Gaia DR3 123456789"]', or None if no match.
    """

    simbad = Simbad()
    simbad.ROW_LIMIT = -1

    # query SIMBAD
    try:
        response = simbad.query_region(targets, radius=radius)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not len(response):
        return RETURNS.NULL

    # convert to DataFrame
    df = response.to_pandas()

    # remove extra whitespace from SIMBAD IDs
    ids = np.asarray([re.sub(r"\s+", " ", str(x)) for x in df["main_id"]])

    # SIMBAD sky positions
    simbad_coords = SkyCoord(ra=df["ra"].to_numpy() * u.deg, dec=df["dec"].to_numpy() * u.deg, frame="icrs")

    # many-to-many sky match
    target_idx, simbad_idx, sep, _ = search_around_sky(targets, simbad_coords, radius)

    # get all matches per target
    matches = [[] for _ in range(len(targets))]
    for t_idx, s_idx in zip(target_idx, simbad_idx):
        matches[t_idx].append(ids[s_idx])

    # remove duplicates while preserving order
    matches = [list(dict.fromkeys(m)) for m in matches]

    # format as strings
    string_results = [", ".join(f"{x}" for x in m) if m else "None" for m in matches]

    # force one data type for fits saving
    max_len = max(len(s) for s in string_results if s)
    results = np.asarray(string_results, dtype=f"<U{max_len}")

    return results
