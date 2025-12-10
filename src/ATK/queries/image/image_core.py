import warnings
from io import BytesIO

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.hdu import ImageHDU
from astropy.wcs import WCS, FITSFixedWarning

from ...structures.definitions import Image
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request

warnings.simplefilter("ignore", category=FITSFixedWarning)


def get_image_data(survey: str, url: str) -> ImageHDU:
    response = send_request(survey, url)
    if response is RETURNS.EXCEPTION:
        return response

    img = fits.open(BytesIO(response.content))[0]
    hdu = ImageHDU(data=img.data, header=img.header)

    return hdu


def get_image(survey: str, band: str, size: int, focus: SkyCoord, hdu: ImageHDU) -> list[Image]:
    wcs = WCS(hdu.header)

    image = Image(survey=survey, band=band, size=size, focus=focus, hdu=hdu, wcs=wcs)

    return [image]
