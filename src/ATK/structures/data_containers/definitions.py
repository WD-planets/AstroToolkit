from dataclasses import dataclass

from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import PrimaryHDU
from astropy.wcs import WCS


@dataclass
class Image:
    survey: str | None = "panstarrs"
    band: str | None = "g"
    hdu: PrimaryHDU | None = None
    wcs: WCS | None = None
    focus: SkyCoord | None = None

    def __repr__(self):
        return f"{self.survey} {self.band}-band Image"

    def __str__(self):
        return self.__repr__()
