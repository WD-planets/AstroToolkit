from dataclasses import dataclass
from pathlib import Path

import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.time import Time
from astropy.wcs import WCS

from ..io.files.writing import write_local
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

    def show(self, show_all_types=False) -> None:
        pprint_structure(self, show_all_types)

    def save(self, path: str | Path = None) -> Path:
        return write_local(self, path)

    def __repr__(self):
        return f"<{self.survey} {self.kind} data>"

    def __str__(self):
        return self.__repr__()

    @property
    def _fname(self):
        if self.source:
            return Path(f"{self.source}_{self.survey}_ATKdata.fits")
        elif self.position:
            return Path(f"{self.position.ra:.3f}_{self.position.dec:.3f}_ATKdata.fits")
        else:
            raise ValueError("No source or position data to be used in generating a file name.")


@dataclass
class Lightcurve:
    survey: str
    band: str
    mjd: numpy.ndarray
    flux: numpy.ndarray | None = None
    flux_err: numpy.ndarray | None = None
    mag: numpy.ndarray | None = None
    mag_err: numpy.ndarray | None = None

    def __post_init__(self):
        # check for a valid input combination
        if (self.flux is None) == (self.mag is None):
            raise ValueError("Lightcurve container cannot hold both 'mag' and 'flux'.")
        if self.flux is not None and self.mag_err is not None:
            raise ValueError("Lightcurve cannot hold invalid combination of 'flux' and 'mag_err'.")
        if self.mag is not None and self.flux_err is not None:
            raise ValueError("Lightcurve cannot hold invalid combination of 'mag' and 'flux_err'.")

        # delete unneeded attributes
        if self.flux is None:
            del self.flux
            del self.flux_err
        if self.mag is None:
            del self.mag
            del self.mag_err

    def __repr__(self):
        return f"<{self.survey} {self.band}-band Lightcurve>"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pandas.DataFrame:
        return struct_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return struct_to_hdu(self)


@dataclass
class Image:
    survey: str | None = None
    band: str | None = None
    hdu: PrimaryHDU | None = None
    wcs: WCS | None = None
    focus: SkyCoord | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band Image>"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pandas.DataFrame:
        return struct_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return struct_to_hdu(self)


@dataclass
class Spectrum:
    survey: str | None = None
    wavelength: numpy.ndarray | None = None
    flux: numpy.ndarray | None = None

    def __repr__(self):
        return f"<{self.survey} Spectrum>"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pandas.DataFrame:
        return struct_to_dataframe(self)

    def to_hdu(self) -> BinTableHDU:
        return struct_to_hdu(self)
