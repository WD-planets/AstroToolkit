from dataclasses import dataclass
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .structures_core import Container
from .Target import Target


@dataclass(repr=False)
class HRD(Container):
    """
    Container for storing the location of sources on the Hertzsprung-Russell Diagram. This object stores both data and relevant metadata.
    """

    #: DOC_OVERRIDE
    survey: str | None = None
    #: Gaia Source ID.
    identifier: int | None = None
    #: DOC_OVERRIDE
    correction: str | None = None
    #: Apparent magnitude band from which to calculate absolute magnitude.
    #:
    #: E.g. ``abs_mag_band = 'Gmag'``
    abs_mag_band: str | None = None
    #: Colour bands.
    #:
    #: E.g. ``colour_bands = 'BPmag-RPmag'``
    colour_bands: str | None = None

    #: Colour value in colour defined by ``colour_bands``.
    colour: numpy.ndarray | None = None
    #: Absolute magnitude value in band given by ``abs_mag_band``, calculated using ``distance``.
    abs_mag: numpy.ndarray | None = None
    #: Distance to source, as derived from parallax.
    distance: Quantity | None = None

    _required = ["survey"]

    _units = {"colour": u.mag, "abs_mag": u.mag, "distance": u.pc}

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def __repr__(self):
        return f"<{self.survey} {self.abs_mag_band} vs {self.colour_bands} HRD>"

    @classmethod
    def from_dataframe(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_dataframe(target, data, **kwargs)

    from_dataframe.__func__.__doc__ = get_docstring("from_dataframe", obj="HRD", args=", ".join(f"``{p}``" for p in _required))

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_table(target, data, **kwargs)

    from_table.__func__.__doc__ = get_docstring("from_table", obj="HRD", args=", ".join(f"``{p}``" for p in _required))
