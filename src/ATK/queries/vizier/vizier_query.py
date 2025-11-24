import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from ...utilities.targeting import check_targeting

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
    except TimeoutError:
        return -1

    if not data:
        return None

    return data[0].to_pandas().sort_values(by=["_r"]).reset_index(drop=True)


def query_by_source(source: int) -> pd.DataFrame | None | int:
    v = Vizier(columns=["**"], column_filters={"Source": f"=={source}"}, row_limit=ROW_LIMIT)

    try:
        data = v.query_constraints(catalog="I/355/gaiadr3", Source=source)
    except TimeoutError:
        return -1

    if not data:
        return None

    return data[0].to_pandas().reset_index(drop=True)


def query(target: SkyCoord | int) -> pd.DataFrame | None | int:
    search_pos, source = check_targeting(target)

    if source:
        return query_by_source(source)
    elif search_pos:
        return query_by_position(search_pos)
    else:
        raise ValueError("query() received no targeting information.")


if __name__ == "__main__":
    target = SkyCoord(ra=141.1853 * u.deg, dec=8.0308 * u.deg, frame="icrs")
    data = query_by_position(target, 300, "I/355/gaiadr3")
