from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.wcs import WCS

from .container_io import container_to_dataframe, container_to_hdu


@dataclass
class Image:
    survey: str | None = "panstarrs"
    band: str | None = "g"
    hdu: PrimaryHDU | None = None
    wcs: WCS | None = None
    focus: SkyCoord | None = None
    test_arr: np.ndarray | None = None

    def __repr__(self):
        return f"{self.survey} {self.band}-band Image"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pd.DataFrame:
        return container_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return container_to_hdu(self)


@dataclass
class Spectrum:
    survey: str | None = None
    wavelength: np.ndarray | None = None
    flux: np.ndarray | None = None

    def __repr__(self):
        return f"{self.survey} Spectrum"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pd.DataFrame:
        return container_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return container_to_hdu(self)
