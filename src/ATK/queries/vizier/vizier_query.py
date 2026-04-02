import warnings

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.exceptions import NoResultsWarning
from astroquery.vizier import Vizier

from ...configuration.survey_config import SURVEY_CONFIG
from ...structures.Record import Record
from ...structures.Target import Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS

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
        data = v.query_region(position, width=radius, catalog=catalogue)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    df = data[0].to_pandas().sort_values(by=["_r"]).reset_index(drop=True)

    return df


def gaia_query_by_source(source: int, kind="data") -> pd.DataFrame | RETURNS:
    """
    Perform a vizier query to Gaia DR3 by source_id
    """

    if kind == "data":
        catalogue = "I/355/gaiadr3"
    elif kind == "lightcurve":
        catalogue = "I/355/epphot"
    else:
        raise ValueError(f"Unexpected Gaia query by source kind '{kind}'.")

    v = Vizier(columns=["**"], column_filters={"Source": f"=={source}"}, row_limit=ROW_LIMIT)

    try:
        data = v.query_constraints(catalog=catalogue, Source=source)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    df = data[0].to_pandas().reset_index(drop=True)

    return df


def query(target: Target, **kwargs) -> pd.DataFrame | RETURNS:
    """
    Perform a Vizier query by source or position (source will be present in kwargs) in the latter case
    """

    aliases = SURVEY_CONFIG._get_aliases()

    survey = kwargs.get("survey")

    # try to get catalogue from alias file
    if survey and survey not in aliases:
        catalogue = survey
    elif survey:
        catalogue = aliases[survey]

    # perform source query if a source was provided
    if target.identifier and catalogue == "I/355/gaiadr3":
        df = gaia_query_by_source(target.identifier)
    else:
        # otherwise perform query by position
        df = query_by_position(target.coords, kwargs["radius"], catalogue)

    if df is RETURNS.NULL or df is RETURNS.EXCEPTION:
        return df

    return [Record(survey=survey, catalogue=catalogue, search_pos=target.coords, data=df, correction=target.correction)]
