from dataclasses import dataclass

import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits import ImageHDU
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs import WCS

from .structures_core import Container


@dataclass(repr=False)
class Image(Container):
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    band: str | None = None

    size: Quantity | None = None
    epoch: Time | None = None
    hdu: ImageHDU | None = None
    wcs: WCS | None = None
    overlay: pandas.DataFrame | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import image_to_hdu

        return image_to_hdu(self)
