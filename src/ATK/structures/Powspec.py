from dataclasses import dataclass
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .structures_core import (Container, DataFrameIOMixin, FITSIOMixin,
                              TableIOMixin, manage_inplace)
from .Target import Target


@dataclass(repr=False)
class Powspec(Container, DataFrameIOMixin, TableIOMixin, FITSIOMixin):
    """
    Container for storing a Lomb-Scargle periodogram. This object stores both data and relevant metadata.
    """

    #: DOC_OVERRIDE
    survey: str | None = None
    #: DOC_OVERRIDE
    band: str | None = None
    #: Survey-specific object ID.
    #:
    #: ``None`` unless present from original :class:`~ATK.Models.Lightcurve`.
    obj_id: str | None = None
    _multiband: bool | None = None
    #: Frequency values.
    frequency: Quantity | None = None
    #: Power at each frequency.
    power: numpy.ndarray | None = None
    #: Frequency of highest power.
    fopt: Quantity | None = None
    #: Period corresponding to frequency of highest power.
    popt: Quantity | None = None

    _data_methods: tuple = ("crop",)

    _required = ["survey", "band"]

    _units = {"fopt": 1 / u.day, "popt": u.day, "frequency": 1 / u.day}

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.power]

        x, ys = crop_nd(x=struct.frequency, ys=ys, lower_lim=min, upper_lim=max)

        struct.frequency = x
        struct.power = ys[0]

        return struct

    crop.__doc__ = get_docstring("crop", x="``frequency``", name="Powspec")
