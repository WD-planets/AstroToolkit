from dataclasses import dataclass, field

import astropy.units as u
import numpy
from astropy.units import Quantity

from ..configuration.base_config import BASE_CONFIG
from .structures_core import Container, manage_inplace

default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")


@dataclass(repr=False)
class SED(Container):
    survey: numpy.ndarray | None = None
    correction: numpy.ndarray | None = None
    band: numpy.ndarray | None = None
    id: numpy.ndarray | None = None
    separation: numpy.ndarray | Quantity | None = None
    wavelength: numpy.ndarray | Quantity | None = None
    flux: numpy.ndarray | Quantity | None = None
    flux_err: numpy.ndarray | Quantity | None = None

    _data_methods: tuple = ("crop", "bin")

    _units = {"flux": u.mJy, "flux_err": u.mJy, "separation": default_scale, "wavelength": u.angstrom}

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def __repr__(self):
        return "<Spectral Energy Distribution>"

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux, struct.flux_err, struct.survey, struct.correction, struct.band, struct.separation]

        x, ys = crop_nd(x=struct.wavelength, ys=ys, lower_lim=min, upper_lim=max)

        flux, flux_err, survey, correction, band, separation = ys

        struct.wavelength = x
        struct.flux = flux
        struct.flux_err = flux_err
        struct.survey = survey
        struct.correction = correction
        struct.band = band
        struct.separation = separation

        return struct
