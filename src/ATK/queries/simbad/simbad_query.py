import re

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad

from ...utilities.defaults import CONNECTION_ERRORS, RETURNS


def get_ids(targets: SkyCoord, radius: float):
    """
    Fetches a list of SIMBAD IDs or None for an astropy SkyCoord (which may contain multiple sets of coordinates)
    """

    simbad = Simbad()
    simbad.ROW_LIMIT = -1

    # send simbad query
    try:
        response = simbad.query_region(targets, radius=radius * u.arcsec)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not response:
        return RETURNS.NULL

    # convert output to pandas dataframe and extract IDs
    df = response.to_pandas()
    ids = [str(id) for id in df["main_id"].tolist()]

    # remove additional whitespace in IDs
    ids = np.asarray([re.sub(r"\s+", " ", s) for s in ids])

    # create SkyCoord from returned SIMBAD positions
    simbad_coords = SkyCoord(ra=df["ra"].to_numpy() * u.deg, dec=df["dec"].to_numpy() * u.deg, frame="icrs")

    # match positions of returned SIMBAD IDs to original targets
    index, separation, _ = simbad_coords.match_to_catalog_sky(targets)

    # keep only matches within requested radius (array of Bools)
    valid = separation <= radius * u.arcsec

    # set up results array
    results = ["None"] * len(targets)

    # iterate through valid IDs, get index of target in original SkyCoord and set this index in results to the returned ID
    for simbad_index in np.where(valid)[0]:
        target_index = index[simbad_index]
        results[target_index] = ids[simbad_index]

    return results
