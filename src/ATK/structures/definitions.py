import copy
from dataclasses import dataclass, field
from pathlib import Path

import astropy.units as u
import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import BinTableHDU, ImageHDU
from astropy.time import Time
from astropy.units import Quantity, Unit
from astropy.wcs import WCS
from bokeh.plotting import figure as Figure

# -------------
# QUERY RESULTS
# -------------

# whether to combine data structures into combined plots
PLOT_METHODS = {"image": "individual", "lightcurve": "combined", "spectrum": "individual", "sed": "individual", "hrd": "combined"}
# whether to split plots by target
SPLIT_BY_TARGET = {"image": True, "lightcurve": True, "spectrum": True, "sed": True, "hrd": False}

# type hint for arrays of astropy Quantities
QuantityArray = numpy.ndarray[Quantity]


def manage_inplace(structure: any, inplace: bool):
    if inplace:
        return structure
    else:
        if hasattr(structure, "figure") and structure.figure:
            structure.figure = None

        return copy.deepcopy(structure)


@dataclass
class Target:
    initial_coords: SkyCoord
    coords: SkyCoord

    identifier: int | None = None
    survey: str | None = None
    correction: str = "none"

    _key: str = field(init=False)
    _aliases: set[str] = field(default_factory=set, init=False)

    def __post_init__(self):
        id_key = f"id:{self.identifier}"
        coord_key = f"coord:{self.initial_coords.ra.deg:.8f},{self.initial_coords.dec.deg:.8f}"

        if self.identifier is not None:
            self._key = id_key
            self._aliases.add(id_key)
        else:
            self._key = coord_key
        self._aliases.add(coord_key)

    @property
    def frame(self):
        return self.coords.frame.name

    @property
    def epoch(self):
        return self.coords.obstime.fits

    @property
    def initial_frame(self):
        return self.initial_coords.frame.name

    @property
    def initial_epoch(self):
        return self.initial_coords.obstime.fits

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

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

        icrs_pos = position.transform_to("icrs")

        return cls(copy.deepcopy(icrs_pos), copy.deepcopy(icrs_pos), None, None, "none")


@dataclass
class BaseQueryResult:
    kind: str | None = None
    survey: str | None = None
    targets: list[Target] | None = field(default_factory=list)
    radius: Quantity | None = None
    exception: bool | None = False

    # maps per-Target key to Target
    _key_map: dict[str, Target] = field(init=False, default_factory=dict)
    # maps per-Target alias to per-Target key
    _alias_map: dict[str, str] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self._build_target_maps()

    def __repr__(self):
        return f"<{self.survey} {self.kind} data>"

    def __str__(self):
        return self.__repr__()

    def _build_target_maps(self):
        self._key_map = {t._key: t for t in self.targets}
        self._alias_map = {}
        for t in self.targets:
            for alias in t._aliases:
                self._alias_map[alias] = t._key

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

        return self

    def save(self, path: str | Path = None) -> Path:
        from ..io.files.writing import write_local

        return write_local(self, path)

    def _fetch_by_key(self, key: str):
        return [d for d in self.data if d._target_key == key]

    def fetch_by_id(self, id: int):
        key = self._alias_map.get(f"id:{id}")
        if key is None:
            return []
        return self._fetch_by_key(key)

    def fetch_by_coord(self, coord: SkyCoord, radius: Quantity | None = 3 * u.arcsec):
        for t in self.targets:
            if coord.separation(t.initial_coords) < radius:
                return self._fetch_by_key(t._key)
        return []

    def fetch_by_target(self, target: Target):
        return self._fetch_by_key(target._key)

    @property
    def _fname(self):
        if self.targets[0].identifier:
            fname = f"{self.targets[0].identifier}"
        else:
            fname = f"{self.position.ra.value:.3f}_{self.position.dec.value:.3f}"
        if self.survey:
            fname += f"_{self.survey}"
        if len(self.targets) > 1:
            fname += "_multi"
        fname += f"_ATK{self.kind}.fits"

        return fname

    def apply(self, method: str, inplace=True, *args, **kwargs):
        struct = manage_inplace(self, inplace)

        for ctnr in struct.data:
            if hasattr(ctnr, method):
                getattr(ctnr, method)(*args, **kwargs)
            else:
                raise ValueError(f"{self.kind} data does not support the method '{method}'.")

        return struct


@dataclass(repr=False)
class QueryResult(BaseQueryResult):
    data: list = field(default_factory=list)


@dataclass(repr=False)
class PlottableQueryResult(BaseQueryResult):
    data: list = field(default_factory=list)
    figure: Figure | None = None

    _stored_plot_params: dict | None = None

    @property
    def _title(self):
        return f"ATK {self.kind.upper()}"

    @property
    def _plot_method(self):
        return PLOT_METHODS[self.kind]

    @property
    def _split_by_target(self):
        return SPLIT_BY_TARGET[self.kind]

    def plot(self, kind: str | None = None, **kwargs: any):
        from ..io.plot_io import plot_data

        self.figure = plot_data(kind, self, **kwargs)
        self._stored_plot_params = kwargs

        return self

    def open(self, fname: Path | str | None = None, **kwargs: any):
        from ..io.plot_io import open as open_html

        open_html(self, fname=fname, **kwargs)

        return self


# ---------------
# DATA CONTAINERS
# ---------------


@dataclass
class BaseContainer:
    _target_key: str | None = None

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

    def __repr__(self):
        return f"<{self.survey} {type(self).__name__}>"

    def __str__(self):
        return self.__repr__()

    def _get_attr_value(self, attr: str) -> numpy.ndarray:
        val = getattr(self, attr)
        if val is None:
            raise ValueError(f"Empty or invalid {type(self).__name__} attribute '{attr}'.")

        if isinstance(val, Quantity):
            return val.value
        return val

    def _get_attr_unit(self, attr: str) -> Unit:
        val = getattr(self, attr)
        if val is None:
            raise ValueError(f"Empty or invalid {type(self).__name__} attribute '{attr}'.")

        if isinstance(val, Quantity):
            return val.unit
        return None

    def _get_cols(self):
        from ..io.structure_io import get_cols

        return get_cols(self)

    def to_dataframe(self) -> pandas.DataFrame:
        from ..io.structure_io import struct_to_dataframe

        return struct_to_dataframe(self)

    @classmethod
    def from_dataframe(cls, data: pandas.DataFrame, **kwargs: any):
        from ..io.structure_io import struct_from_dataframe

        return struct_from_dataframe(cls, data, **kwargs)

    def to_hdu(self) -> BinTableHDU:
        from ..io.structure_io import struct_to_hdu

        return struct_to_hdu(self)


@dataclass(repr=False)
class VizierEntry(BaseContainer):
    survey: str | None = None
    catalogue: str | None = None
    correction: str | None = None
    search_pos: str | None = None
    separation: str | None = None

    data: pandas.DataFrame | None = None

    def __repr__(self):
        return f"<{self.survey} ({self.catalogue}) Vizier Data>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)


@dataclass(repr=False)
class Lightcurve(BaseContainer):
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    separation: Quantity | None = None
    band: str | None = None

    obj_id: str | None = None
    mjd: numpy.ndarray | None = None
    flux: numpy.ndarray | QuantityArray | None = None
    flux_err: numpy.ndarray | QuantityArray | None = None
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

        # ensure valid combination of flux/flux_err/mag/mag_err
        if self.flux is None:
            for f in ("flux", "flux_err"):
                self.__dict__.pop(f, None)
        if self.mag is None:
            for f in ("mag", "mag_err"):
                self.__dict__.pop(f, None)

    @property
    def brightness(self):
        return getattr(self, self.brightness_type)

    @property
    def brightness_err(self):
        return getattr(self, f"{self.brightness_type}_err")

    def set_brightness(self, val: numpy.ndarray):
        setattr(self, self.brightness_type, val)

    def set_brightness_err(self, val: numpy.ndarray):
        setattr(self, f"{self.brightness_type}_err", val)

    @property
    def brightness_type(self):
        if self.mag is not None:
            return "mag"
        elif self.flux is not None:
            return "flux"
        else:
            # shouldn't happen due to __post_init__
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.lightcurve.binning import bin_nd

        struct = manage_inplace(self, inplace)

        x, ys, errs = bin_nd(x=struct.mjd, ys=[struct.brightness, struct.ra, struct.dec], errs=[struct.brightness_err], bins=bins, size=size)

        brightness, ra, dec = ys
        brightness_err = errs[0]

        struct.mjd = x
        struct.set_brightness(brightness)
        struct.set_brightness_err(brightness_err)
        struct.ra = ra
        struct.dec = dec

        return struct

    def crop(self, min: float, max: float, inplace=True):
        from .methods.lightcurve.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        x, ys = crop_nd(x=struct.mjd, ys=[struct.brightness, struct.brightness_err, struct.ra, struct.dec], lower_lim=min, upper_lim=max)

        brightness, brightness_err, ra, dec = ys

        struct.mjd = x
        struct.set_brightness(brightness)
        struct.set_brightness_err(brightness_err)
        struct.ra = ra
        struct.dec = dec

        return struct


@dataclass(repr=False)
class Image(BaseContainer):
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


@dataclass(repr=False)
class Spectrum(BaseContainer):
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None

    separation: Quantity | None = None
    exposure: Quantity | None = None
    wavelength: numpy.ndarray | QuantityArray | None = None
    flux: numpy.ndarray | QuantityArray | None = None


@dataclass(repr=False)
class SED(BaseContainer):
    survey: numpy.ndarray | None = None
    correction: numpy.ndarray | None = None
    band: numpy.ndarray | None = None

    separation: numpy.ndarray | QuantityArray | None = None
    wavelength: numpy.ndarray | QuantityArray | None = None
    flux: numpy.ndarray | QuantityArray | None = None
    flux_err: numpy.ndarray | QuantityArray | None = None

    def __repr__(self):
        return "<Spectral Energy Distribution>"


@dataclass(repr=False)
class HRD(BaseContainer):
    survey: str | None = None
    identifier: int | None = None
    correction: str | None = None
    abs_mag_band: str | None = None
    colour_bands: str | None = None

    colour: numpy.ndarray | None = None
    abs_mag: numpy.ndarray | None = None
    distance: Quantity | None = None

    def __repr__(self):
        return f"<{self.survey} {self.abs_mag_band} vs {self.colour_bands} HRD>"
