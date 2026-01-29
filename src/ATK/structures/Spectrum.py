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
    velocity: numpy.ndarray | QuantityArray | None = None
    flux: numpy.ndarray | QuantityArray | None = None

    _required: tuple[str] = "flux"

    def __post_init__(self):
        # check for a valid input combination
        if (self.wavelength is None) == (self.velocity is None):
            raise ValueError("Spectrum container must hold one of 'wavelength' and 'velocity'.")

        # ensure valid combination of wavelength/velocity
        if self.velocity is None:
            self.__dict__.pop("velocity", None)
        if self.wavelength is None:
            self.__dict__.pop("wavelength", None)

    @property
    def x_type(self):
        if self.wavelength is not None:
            return "wavelength"
        elif self.velocity is not None:
            return "velocity"
        else:
            raise ValueError("Spectrum container must hold one of 'wavelength' and 'velocity'.")

    @property
    def x_arr(self):
        return getattr(self, self.x_type)

    def set_x(self, val: numpy.ndarray):
        setattr(self, self.x_type, val)

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys = crop_nd(x=struct.x_arr, ys=ys, lower_lim=min, upper_lim=max)

        struct.set_x(x)
        struct.flux = ys[0]

        return struct

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.binning import bin_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys, _ = bin_nd(x=struct.x_arr, ys=ys, errs=[], bins=bins, size=size)

        struct.set_x(x)
        struct.flux = ys[0]

        return struct

    def rv_fit(self, feature_wavelengths: list[float], feature_widths: list[float]):
        from .methods.spectrum.fitting import do_fitting

        do_fitting(self, feature_wavelengths, feature_widths)

    def vspec(self, wav_ref: float | Quantity, inplace: bool = True):
        from .methods.spectrum.fitting import get_velocities

        struct = manage_inplace(self, inplace)

        struct.velocity = get_velocities(struct.wavelength, wav_ref)
        struct.wavelength = None

        return struct
