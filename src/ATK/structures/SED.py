from dataclasses import dataclass, field
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..configuration.base_config import BASE_CONFIG
from ..utilities.docstrings import get_docstring
from .Base import (Container, DataFrameIOMixin, FITSIOMixin, TableIOMixin,
                   manage_inplace)
from .Target import Target

default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")

PLOT_PARAMS = {
    "overlay": [
        ":class:`~ATK.Models.Spectrum`, optional",
        """
        Overlays a spectrum.

        By default, no spectrum is overlayed.
        """,
    ]
}


@dataclass(repr=False)
class SED(Container, DataFrameIOMixin, TableIOMixin, FITSIOMixin):
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

    _data_methods: tuple = ("crop",)

    _required = ["survey"]

    _units = {"flux": u.mJy, "flux_err": u.mJy, "separation": default_scale, "wavelength": u.angstrom}

    _plot_params = PLOT_PARAMS

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def __repr__(self):
        return "<Spectral Energy Distribution>"

    def crop(self, cmin: float | None = None, cmax: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux, struct.flux_err, struct.survey, struct.correction, struct.band, struct.separation]

        x, ys = crop_nd(x=struct.wavelength, ys=ys, lower_lim=cmin, upper_lim=cmax)

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
