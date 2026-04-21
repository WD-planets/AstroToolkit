from dataclasses import dataclass, field
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..configuration.base_config import BASE_CONFIG
from ..utilities.docstrings import get_docstring
from .structures_core import Container, manage_inplace
from .Target import Target

default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")


@dataclass(repr=False)
class SED(Container):
    """
    Container for storing spectral energy distribution data. This object stores both data and relevant metadata.
    """

    #: DOC_OVERRIDE
    survey: numpy.ndarray | None = None
    #: DOC_OVERRIDE
    correction: numpy.ndarray | None = None
    #: Photometric bands.
    band: numpy.ndarray | None = None
    #: Per-survey object IDs.
    id: numpy.ndarray | None = None
    #: Separation between position of the search and the returned photometric detections.
    separation: Quantity | None = None
    #: Wavelength values.
    wavelength: Quantity | None = None
    #: Flux values.
    flux: Quantity | None = None
    #: Flux error values.
    flux_err: Quantity | None = None

    _data_methods: tuple = ("crop", "bin")

    _required = ["survey"]

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

    crop.__doc__ = get_docstring("bin", x="``wavelength``", name="SED")

    @classmethod
    def from_dataframe(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_dataframe(target, data, **kwargs)

    from_dataframe.__func__.__doc__ = get_docstring("from_dataframe", obj="SED", args=", ".join(f"``{p}``" for p in _required))

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_table(target, data, **kwargs)

    from_table.__func__.__doc__ = get_docstring("from_table", obj="SED", args=", ".join(f"``{p}``" for p in _required))
