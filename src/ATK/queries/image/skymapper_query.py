from io import BytesIO

import pandas as pd

from ...structures.definitions import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request
from .image_core import get_image_data, mjd_to_epoch


def check_inputs(band: str, size: int):
    if size > 600:
        raise ValueError("Size too large. Maximum supported by panstarrs is 1500 arcsec.")
    if band not in ["g", "r", "i", "z", "u", "v"]:
        raise ValueError("Invalid panstarrs bands. Supported bands are ['g', 'r', 'i', 'z', 'u', 'v'].")


# size: int, band: str, overlays: list | dict, search_pos: SkyCoord = None, **kwargs: any
def query(target: Target, **kwargs: any):
    band, size = kwargs["band"], kwargs["size"]

    check_inputs(band, size)

    # needed in degrees
    url_size = size / 3600

    # fetch table
    url = f"https://api.skymapper.nci.org.au/public/siap/dr4/query?POS={target.coords.ra.value},{target.coords.dec.value}&SIZE={url_size}&BAND={band}&FORMAT=image/fits&VERB=3&INTERSECT=covers&RESPONSEFORMAT=CSV"

    response = send_request("skymapper", url)

    if response is RETURNS.EXCEPTION:
        return response

    if "text/csv" not in response.headers.get("Content-Type"):
        return RETURNS.NULL

    table = pd.read_csv(BytesIO(response.content))

    if not len(table):
        return RETURNS.NULL

    # get only URLs matching the requested band
    table = table[table["band"] == band]
    if not len(table):
        return RETURNS.NULL

    # sort to get image of maximum exposure and minimum air mass
    table = table.sort_values(["exptime", "airmass"], ascending=[False, True])

    # get url from table
    image_url = table["get_fits"][0]

    image_hdu = get_image_data(image_url, "panstarrs", band, size, mjd_to_epoch, "MJD-OBS")

    return image_hdu
