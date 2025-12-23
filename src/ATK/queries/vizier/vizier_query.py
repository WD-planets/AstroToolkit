import warnings
from urllib.error import HTTPError

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.exceptions import NoResultsWarning
from astroquery.vizier import Vizier
from requests.exceptions import ConnectionError, ConnectTimeout

from ...configuration.alias_config import ALIAS_CONFIG
from ...structures.definitions import Target
from ...utilities.defaults import RETURNS

warnings.simplefilter("ignore", category=NoResultsWarning)

# ensure all rows are returned
ROW_LIMIT = -1
Vizier.ROW_LIMIT = -1


def query_by_position(position: SkyCoord, radius: float, catalogue: str) -> pd.DataFrame | None | int:
    """
    Returns a DataFrame of Vizier catalogue data within a given radius of a given position, sorted by distance to the target
    """

    v = Vizier(columns=["**"], row_limit=ROW_LIMIT)
    try:
        data = v.query_region(position, width=radius * u.arcsec, catalog=catalogue)
    except (TimeoutError, ConnectionError, ConnectTimeout, HTTPError):
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    return data[0].to_pandas().sort_values(by=["_r"]).reset_index(drop=True)


def gaia_query_by_source(source: int) -> pd.DataFrame | None | int:
    """
    Perform a vizier query to Gaia DR3 by source_id
    """

    v = Vizier(columns=["**"], column_filters={"Source": f"=={source}"}, row_limit=ROW_LIMIT)

    try:
        data = v.query_constraints(catalog="I/355/gaiadr3", Source=source)
    except (TimeoutError, ConnectionError, ConnectTimeout):
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    return data[0].to_pandas().reset_index(drop=True)


def query(target: Target, **kwargs) -> pd.DataFrame | RETURNS:
    """
    Perform a Vizier query by source or position (source will be present in kwargs) in the latter case
    """

    aliases = ALIAS_CONFIG.as_dict()["vizier_aliases"]

    # survey = catalogue alias (here for parity with other query commands), catalogue = actual vizier catalogue ID
    survey, catalogue = kwargs.get("survey"), kwargs.get("catalogue")

    # ensure exactly one of 'survey', 'catalogue' provided
    if survey is None == catalogue is None:
        raise ValueError("Either 'survey' or 'catalogue' required for Vizier queries.")

    # try to get catalogue from alias file
    if survey and survey not in aliases:
        raise ValueError(f"Survey '{survey}' not found in ATK alias file.")
    elif survey:
        catalogue = aliases[survey]

    # perform source query if a source was provided
    if target.identifier and catalogue == "I/355/gaiadr3":
        return gaia_query_by_source(target.identifier)

    # otherwise perform query by position
    return query_by_position(target.coords, kwargs["radius"], catalogue)
