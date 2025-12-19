from urllib.error import HTTPError

import astropy.units as u
from astropy.wcs import WCS
from astroquery.sdss import SDSS
from requests.exceptions import ConnectionError, ConnectTimeout

from ...structures.definitions import Image, Target
from ...utilities.defaults import RETURNS
from .image_core import get_image_skycoord, iso_to_epoch, reproject_hdu


def query(target: Target, **kwargs) -> Image:
    band, size = kwargs["band"], kwargs["size"]

    # get image list
    try:
        imgs = SDSS.get_images(coordinates=target.coords, band=band, radius=size * u.arcsec)
    except (TimeoutError, ConnectionError, ConnectTimeout, HTTPError):
        return RETURNS.EXCEPTION

    if not imgs:
        return RETURNS.NULL

    # get first image
    hdu = imgs[0][0]

    # reproject image to north_up
    proj_hdu = reproject_hdu(hdu, target.coords, size)

    # get information for Image
    image_focus = get_image_skycoord(proj_hdu, iso_to_epoch, "DATE-OBS")
    wcs_out = WCS(proj_hdu.header)

    image = Image("sdss", band, size, proj_hdu, wcs_out, image_focus, image_focus.obstime)

    return [image]
