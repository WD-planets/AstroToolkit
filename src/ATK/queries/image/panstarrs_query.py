from io import BytesIO

import astropy.units as u
from astropy.table import Table

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from ...utilities.misc import suppress_stdout
from ...utilities.requests import send_request
from .image_core import get_image_data, mjd_to_epoch


def check_inputs(band: str, size: int):
    if size > 1500 * u.arcsec:
        raise ValueError("Size too large. Maximum supported by panstarrs is 1500 arcsec.")
    bands = ["g", "r", "i", "z", "y"]
    if band not in bands:
        raise ValueError(f"Invalid panstarrs bands. Supported bands are: {', '.join(bands)}.")


# size: int, band: str, overlays: list | dict, search_pos: SkyCoord = None, **kwargs: any
def query(target: Target, **kwargs: any):
    """
    Performs a PanStarrs image query
    """

    band, size = kwargs["band"], kwargs["size"]

    check_inputs(band, size)

    # 0.25 arcsec per pixel
    url_size = int(size.to(u.arcsec).value * 4)

    # fetch table
    url = f"https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={target.coords.ra.value}&dec={target.coords.dec.value}&band={band}"
    with suppress_stdout():
        response = send_request("panstarrs", url)
    if response is RETURNS.EXCEPTION:
        return response

    table = Table.read(BytesIO(response.content), format="ascii")
    if not len(table):
        return RETURNS.NULL

    # get only table rows matching the requested band
    table = table[table["filter"] == band]
    if not len(table):
        return RETURNS.NULL

    # get url from table
    sub_url = f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={target.coords.ra.value}&dec={target.coords.dec.value}&size={url_size}&format=fits&red="
    fname = table["filename"][0]
    main_url = f"{sub_url}{fname}"

    image = get_image_data(main_url, "panstarrs", band, size, mjd_to_epoch, "MJD-OBS")

    return image
