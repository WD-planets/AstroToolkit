from dataclasses import dataclass
from pathlib import Path

import numpy
import pandas
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.time import Time
from astropy.wcs import WCS

from ..io.files.writing import write_structure
from ..io.struct_stdout import pprint_structure
from .structure_io import struct_to_dataframe, struct_to_hdu


@dataclass
class QueryResult:
    kind: str | None = None
    survey: str | None = None
    radius: float | None = None
    source: int | None = None
    position: SkyCoord | None = None
    epoch: Time | None = None
    frame: str | None = None
    correction: str | None = None
    exception: bool | None = False

    data: pandas.DataFrame | list | None = None

    def show(self, show_all_types=False):
        pprint_structure(self, show_all_types)

    def save(self, path: str | Path = None):
        write_structure(self, path)

    def __repr__(self):
        return f"{self.survey} {self.kind} data"

    def __str__(self):
        return self.__repr__()

    @property
    def _fname(self):
        if self.source:
            return f"{self.source}_{self.survey}_data.fits"
        elif self.position:
            return f"{self.position.ra:.3f}_{self.position.dec:.3f}_{self.kind}.fits"


@dataclass
class Image:
    survey: str | None = "panstarrs"
    band: str | None = "g"
    hdu: PrimaryHDU | None = None
    wcs: WCS | None = None
    focus: SkyCoord | None = None
    test_arr: numpy.ndarray | None = None

    def __repr__(self):
        return f"{self.survey} {self.band}-band Image"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pd.DataFrame:
        return struct_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return struct_to_hdu(self)


@dataclass
class Spectrum:
    survey: str | None = None
    wavelength: numpy.ndarray | None = None
    flux: numpy.ndarray | None = None

    def __repr__(self):
        return f"{self.survey} Spectrum"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pd.DataFrame:
        return struct_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return struct_to_hdu(self)
