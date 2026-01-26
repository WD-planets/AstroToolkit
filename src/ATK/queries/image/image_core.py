import bz2
import gzip
import warnings
from datetime import datetime, timedelta
from io import BytesIO
from types import FunctionType

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits.hdu import ImageHDU
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from reproject import reproject_interp

from ...structures.Image import Image
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request

# ignore fits warning
warnings.filterwarnings("ignore", category=FITSFixedWarning)


def mjd_to_epoch(hdr: Header, key: str) -> Time:
    """
    Convert an mjd in a given header key to an astropy Time
    """

    mjd = hdr.get(key)
    if not mjd:
        raise ValueError(f"Epoch Key '{key}' not found in image header.")

    return Time(mjd, format="mjd")


def iso_to_epoch(hdr: Header, key: str) -> Time:
    """
    Convert an ISO time in a given header key to an astropy Time
    """

    iso = hdr.get(key)
    if not iso:
        raise ValueError(f"Epoch Key '{key}' not found in image header.")

    return Time(iso, format="iso")


def fits_to_epoch(hdr: Header, key: str) -> Time:
    """
    Convert a fits time in a given header key to an astropy Time
    """

    fits_time = hdr.get(key)
    if not fits_time:
        raise ValueError(f"Epoch Key '{key}' not found in image header.")

    # sometimes fits headers from DSS contain invalid minutes (i.e. minutes exceed 60)
    date, time = fits_time.split("T")

    # convert all components of time to integers (at most loses seconds of precision)
    year, month, day = map(int, map(float, date.split("-")))
    hour, minute, second = map(int, map(float, time.split(":")))

    # fix time
    base = datetime(year, month, day, hour, 0, 0)
    delta = timedelta(minutes=minute, seconds=second)
    fixed_time = base + delta

    return Time(fixed_time)


def get_image_skycoord(hdu: ImageHDU, epoch_fetcher: FunctionType | None = None, epoch_key: str | None = None) -> SkyCoord:
    """
    Returns the central SkyCoord of an ImageHDU, optionally fetching the epoch using a function epoch_fetcher.
    This function may take a key 'epoch_key', e.g. MJD-OBS, in which case this will be fetched from the header before transformation.
    If no key is provided, epoch_fetcher is assumed to only take the HDU as an argument
    """

    ny, nx = hdu.data.shape

    wcs = WCS(hdu.header)

    # central pixel
    x_center = (nx - 1) / 2
    y_center = (ny - 1) / 2

    # get coordinates of central pixel
    centre = wcs.pixel_to_world(x_center, y_center)

    # if no epoch was given, just return central coordinates
    if not epoch_fetcher:
        return centre

    # standard epoch extraction functions
    if epoch_key:
        obstime = epoch_fetcher(hdu.header, epoch_key)

    # non-standard (uses hard-coded keys)
    else:
        obstime = epoch_fetcher(hdu.header)

    centre_skycoord = SkyCoord(centre.ra, centre.dec, frame=centre.frame, obstime=obstime).transform_to("icrs")

    return centre_skycoord


def get_image_size(hdu: ImageHDU) -> tuple[float]:
    """
    Returns the actual size of an image as tuple(x,y)
    """

    wcs = WCS(hdu.header)

    # Image shape
    ny, nx = hdu.data.shape

    # Pixel scales in degrees/pixel
    scales = proj_plane_pixel_scales(wcs)

    # Compute total angular size
    size_x = round(nx * scales[0] * 3600, 1)
    size_y = round(ny * scales[1] * 3600, 1)

    size = (size_x * u.arcsec, size_y * u.arcsec)

    return size


def get_image_data(url: str, survey: str, band: str, size: int, epoch_fetcher: FunctionType, epoch_key: str) -> list[ImageHDU]:
    """
    Returns a list (for parity) containing an ATK image with all returned data. May also return one of RETURNS
    """

    response = send_request(survey, url)
    if response is RETURNS.EXCEPTION:
        return response

    # decompress fits files if needed
    if response.content.startswith(b"\x1f\x8b"):
        # gzip - response starts with bytes '1F 8B'
        data = gzip.decompress(response.content)
    elif response.content.startswith(b"BZh"):
        # bzip2 - response starts with BZh
        data = bz2.decompress(response.content)
    else:
        data = response.content

    img = fits.open(BytesIO(data))[0]

    hdu = ImageHDU(data=img.data, header=img.header)
    img_centre = get_image_skycoord(img, epoch_fetcher, epoch_key)
    wcs = WCS(hdu.header)

    image = Image(survey=survey, band=band, size=size, search_pos=img_centre, hdu=hdu, wcs=wcs, epoch=epoch_fetcher(hdu.header, epoch_key))

    return [image]


def reproject_hdu(hdu: ImageHDU, cutout_centre: SkyCoord, cutout_size: float) -> ImageHDU:
    """
    Reprojects a fits ImageHDU to north-up, and returns a hdu containing a cutout of size cutout_size (arcsec) around cutout_centre
    """

    data = hdu.data.astype(float)
    wcs = WCS(hdu.header)
    header = hdu.header.copy()

    # get pixel sizes in x and y
    pixel_scale_x, pixel_scale_y = [s.to(u.arcsec).value for s in wcs.proj_plane_pixel_scales()]

    # calculate output array sizes in pixels
    n_pix_x = int(cutout_size / pixel_scale_x)
    n_pix_y = int(cutout_size / pixel_scale_y)

    # ensure odd-shaped arrays to centre target
    if n_pix_x % 2 == 0:
        n_pix_x += 1
    if n_pix_y % 2 == 0:
        n_pix_y += 1

    # create output WCS, squared and north-up, centred on cutout_centre
    wcs_out = WCS(naxis=2)
    wcs_out.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_out.wcs.crval = [cutout_centre.ra.deg, cutout_centre.dec.deg]
    wcs_out.wcs.crpix = [(n_pix_x + 1) / 2, (n_pix_y + 1) / 2]
    wcs_out.wcs.cdelt = [-pixel_scale_x / 3600.0, pixel_scale_y / 3600.0]

    # reproject using nearest-neighbour matching
    data_proj, _ = reproject_interp((data, wcs), wcs_out, shape_out=(n_pix_y, n_pix_x), order=0)

    # ensure array is writeable
    data_proj = np.array(data_proj, copy=True)

    # merge original header with new WCS
    new_header = header.copy()
    for key, val in wcs_out.to_header().items():
        new_header[key] = val

    return ImageHDU(data=data_proj, header=new_header)
