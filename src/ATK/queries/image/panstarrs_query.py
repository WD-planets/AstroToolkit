from io import BytesIO

from astropy.table import Table

from ...structures.definitions import Target
from ...utilities.defaults import RETURNS
from ...utilities.misc import suppress_stdout
from ...utilities.requests import send_request
from .image_core import get_image_data, mjd_to_epoch


def check_inputs(band: str, size: int):
    if size > 1500:
        raise ValueError("Size too large. Maximum supported by panstarrs is 1500 arcsec.")
    if band not in ["g", "r", "i", "z", "y"]:
        raise ValueError("Invalid panstarrs bands. Supported bands are ['g', 'r', 'i', 'z', 'y'].")


# size: int, band: str, overlays: list | dict, search_pos: SkyCoord = None, **kwargs: any
def query(target: Target, **kwargs: any):
    band, size = kwargs["band"], kwargs["size"]

    check_inputs(band, size)

    # 0.25 arcsec per pixel
    url_size = size * 4

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

    image_hdu = get_image_data(main_url, "panstarrs", band, size, mjd_to_epoch, "MJD-OBS")

    return image_hdu
