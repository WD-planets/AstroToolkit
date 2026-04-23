from dataclasses import dataclass
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .Base import Container, DataFrameIOMixin, FITSIOMixin, TableIOMixin
from .Target import Target

PLOT_PARAMS = {
    "split": [
        "bool, optional",
        """
        If True, multiple HRDs in a :class:`~ATK.Models.DataSet` are plotted independently. Otherwise, combine into a single HRD plot."

        Default is ``True``.
        """,
    ],
    "background": [
        "float, optional",
        """
        Sets the fraction of the background sample that is rendered, from ``0.0`` for no background to ``1.0`` for the full background sample.

        Default is ``1.0``.
        """,
    ],
}


@dataclass(repr=False)
class HRD(Container, DataFrameIOMixin, TableIOMixin, FITSIOMixin):
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

    _plot_params = PLOT_PARAMS

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def __repr__(self):
        return f"<{self.survey} {self.abs_mag_band} vs {self.colour_bands} HRD>"
