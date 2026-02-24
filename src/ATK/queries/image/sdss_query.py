import astropy.units as u
from astropy.wcs import WCS
from astroquery.sdss import SDSS

from ...structures.Image import Image
from ...structures.Target import Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from .image_core import get_image_skycoord, iso_to_epoch, reproject_hdu


def check_inputs(band: str, size: int):
    if size > 500 * u.arcsec:
        raise ValueError("Size too large. Maximum supported by panstarrs is 1500 arcsec.")
    if band not in ["u", "g", "r", "i", "z"]:
        raise ValueError("Invalid panstarrs bands. Supported bands are ['u', 'g', 'r', 'i', 'z'].")


def query(target: Target, **kwargs) -> Image:
    """
    Performs an SDSS image query
    """

    band, size = kwargs["band"], kwargs["size"]

    check_inputs(band, size)

    # get image list
    try:
        imgs = SDSS.get_images(coordinates=target.coords, band=band, radius=size)
    except CONNECTION_ERRORS:
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
