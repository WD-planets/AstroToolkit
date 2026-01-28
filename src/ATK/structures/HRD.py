from dataclasses import dataclass

import numpy
from astropy.units import Quantity

from .structures_core import Container


@dataclass(repr=False)
class HRD(Container):
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
