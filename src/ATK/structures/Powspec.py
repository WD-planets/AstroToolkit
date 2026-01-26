from dataclasses import dataclass

import numpy
from astropy.units import Quantity

from .structures_core import BaseContainer
from .structures_core import QuantityArray, manage_inplace


@dataclass(repr=False)
class Powspec(BaseContainer):
    survey: str | None = None
    band: str | None = None

    obj_id: str | None = None
    frequency: QuantityArray | None = None
    power: numpy.ndarray | None = None
    fopt: Quantity | None = None
    popt: Quantity | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def crop(self, min: float, max: float, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.power]

        x, ys = crop_nd(x=struct.frequency, ys=ys, lower_lim=min, upper_lim=max)

        struct.frequency = x
        struct.power = ys[0]

        return struct
