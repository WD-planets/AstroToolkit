from dataclasses import dataclass, field
from pathlib import Path

import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import BinTableHDU, ImageHDU
from astropy.time import Time
from astropy.wcs import WCS
from bokeh.plotting import figure as Figure

# -------------
# QUERY RESULTS
# -------------

PLOT_METHODS = {"image": "individual", "lightcurve": "combined", "spectrum": "individual", "sed": "individual"}


@dataclass
class Target:
    coords: SkyCoord
    identifier: int | None = None
    survey: str | None = None
    correction: str = "none"

    def show(self, show_all_types=False) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types)

    @classmethod
    def from_id(cls, id: int, survey="gaia"):
        if survey == "gaia":
            from ..utilities.coordinates import get_gaia_target

            return get_gaia_target(id)
        else:
            raise NotImplementedError("Other astronometric surveys will be added at a later date.")

    @classmethod
    def from_pos(cls, position: SkyCoord):
        # if no epoch was set, assume J2000
        if not position.obstime:
            j2000 = Time("2000-01-01T00:00:00.000", format="fits")

            position = SkyCoord(position.data, frame=position.frame, obstime=j2000)

        return cls(position.transform_to("icrs"), None, None, "none")


@dataclass
class BaseQueryResult:
    kind: str | None = None
    survey: str | None = None
    radius: float | None = None
    position: SkyCoord | None = None
    identifier: int | None = None
    epoch: Time | None = None
    frame: str | None = None
    correction: str | None = None
    exception: bool | None = False

    def show(self, show_all_types=False) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types)

    def save(self, path: str | Path = None) -> Path:
        from ..io.files.writing import write_local

        return write_local(self, path)

    def __repr__(self):
        return f"<{self.survey} {self.kind} data>"

    def __str__(self):
        return self.__repr__()

    @property
    def _fname(self):
        if self.identifier:
            return Path(f"{self.identifier}_{self.survey}_ATKdata.fits")
        elif self.position:
            return Path(f"{self.position.ra.value:.3f}_{self.position.dec.value:.3f}_{self.survey}_ATKdata.fits")
        else:
            raise ValueError("No source or position data to be used in generating a file name.")


@dataclass(repr=False)
class QueryResult(BaseQueryResult):
    data: list = field(default_factory=list)


@dataclass(repr=False)
class PlottableQueryResult(BaseQueryResult):
    data: list = field(default_factory=list)
    figure: Figure | None = None

    @property
    def _plot_method(self):
        return PLOT_METHODS[self.kind]

    def plot(self, kind: str | None = None, **kwargs: any):
        from .plot_io import plot_data

        self.figure = plot_data(kind, self, **kwargs)

    def open(self, fname: Path | str | None = None):
        from .plot_io import open

        open(self, fname=fname)


# ---------------
# DATA CONTAINERS
# ---------------


@dataclass
class BaseContainer:
    def show(self, show_all_types=False) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types)

    def __repr__(self):
        return f"<{self.survey} {type(self).__name__}>"

    def __str__(self):
        return self.__repr__()

    def to_dataframe(self) -> pandas.DataFrame:
        from .structure_io import struct_to_dataframe

        return struct_to_dataframe(self)

    @classmethod
    def from_dataframe(cls, data: pandas.DataFrame, **kwargs: any):
        from .structure_io import struct_from_dataframe

        print(kwargs)

        return struct_from_dataframe(cls, data, **kwargs)

    def to_hdu(self) -> BinTableHDU:
        from .structure_io import struct_to_hdu

        return struct_to_hdu(self)


@dataclass(repr=False)
class Lightcurve(BaseContainer):
    survey: str
    band: str
    mjd: numpy.ndarray
    flux: numpy.ndarray | None = None
    flux_err: numpy.ndarray | None = None
    mag: numpy.ndarray | None = None
    mag_err: numpy.ndarray | None = None
    ra: numpy.ndarray | None = None
    dec: numpy.ndarray | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def __post_init__(self):
        # check for a valid input combination
        if (self.flux is None) == (self.mag is None):
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")
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

    @property
    def brightness_type(self):
        if hasattr(self, "flux"):
            return "flux"
        elif hasattr(self, "mag"):
            return "mag"
        else:
            # shouldn't happen due to __post_init__
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")


@dataclass(repr=False)
class Image(BaseContainer):
    survey: str | None = None
    band: str | None = None
    size: int | None = None
    hdu: ImageHDU | None = None
    wcs: WCS | None = None
    focus: SkyCoord | None = None
    epoch: Time | None = None
    overlay: pandas.DataFrame | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from .structure_io import image_to_hdu

        return image_to_hdu(self)


@dataclass(repr=False)
class Spectrum(BaseContainer):
    survey: str | None = None
    position: SkyCoord | None = None
    separation: float | None = None
    exposure: float | None = None
    wavelength: numpy.ndarray | None = None
    flux: numpy.ndarray | None = None


@dataclass(repr=False)
class SED(BaseContainer):
    survey: numpy.ndarray | None = None
    band: numpy.ndarray | None = None
    separation: numpy.ndarray | None = None
    wavelength: numpy.ndarray | None = None
    flux: numpy.ndarray | None = None
    flux_err: numpy.ndarray | None = None

    def __repr__(self):
        return "<Spectral Energy Distribution>"
