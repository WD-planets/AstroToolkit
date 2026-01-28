from dataclasses import dataclass

import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from .structures_core import Container, QuantityArray, manage_inplace


@dataclass(repr=False)
class Spectrum(Container):
    # --- metadata ---
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    separation: Quantity | None = None
    exposure: Quantity | None = None

    # --- data ---
    wavelength: numpy.ndarray | QuantityArray | None = None
    flux: numpy.ndarray | QuantityArray | None = None

    _required: tuple[str] = ("wavelength", "flux")

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
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

    def rv_fit(self, feature_wavelengths: list[float], feature_widths: list[float]):
        from .methods.spectrum.fitting_2 import do_fitting

        do_fitting(self, feature_wavelengths, feature_widths)
