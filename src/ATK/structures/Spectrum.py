from dataclasses import dataclass

import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from .structures_core import BaseContainer
from .structures_core import QuantityArray, manage_inplace


@dataclass(repr=False)
class Spectrum(BaseContainer):
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None

    separation: Quantity | None = None
    exposure: Quantity | None = None
    wavelength: numpy.ndarray | QuantityArray | None = None
    flux: numpy.ndarray | QuantityArray | None = None

    def crop(self, min: float, max: float, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys = crop_nd(x=struct.wavelength, ys=ys, lower_lim=min, upper_lim=max)

        struct.wavelength = x
        struct.flux = ys[0]

        return struct

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.binning import bin_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys, _ = bin_nd(x=struct.wavelength, ys=ys, errs=[], bins=bins, size=size)

        struct.wavelength = x
        struct.flux = ys[0]

        return struct
