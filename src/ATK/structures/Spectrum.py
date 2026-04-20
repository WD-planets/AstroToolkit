from dataclasses import dataclass, field

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from .methods.spectrum.fitting import do_fitting
from .methods.spectrum.radial_velocities import get_rvs
from .structures_core import Container, manage_inplace


@dataclass(repr=False)
class Spectrum(Container):
    # --- metadata ---
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    separation: Quantity | None = None
    exposure: Quantity | None = None
    wav_ref: Quantity | None = None

    # --- data ---
    wavelength: numpy.ndarray | Quantity | None = None
    velocity: numpy.ndarray | Quantity | None = None
    flux: numpy.ndarray | Quantity | None = None

    _required = ["survey"]

    _data_methods: tuple = ("crop", "bin", "vspec")
    _plot_methods: dict = field(default_factory=lambda: {"fit": do_fitting, "rv_fit": get_rvs})

    _units = {
        "exposure": u.s,
        "wav_ref": u.angstrom,
        "wavelength": u.angstrom,
        "velocity": u.km / u.s,
        "flux": u.Unit(1e-17) * u.erg / (u.s * u.cm**2 * u.AA),
    }

    def __post_init__(self):
        # check for a valid input combination
        if (self.wavelength is None) == (self.velocity is None):
            raise ValueError("Spectrum container must hold one of 'wavelength' and 'velocity'.")

        # ensure valid combination of wavelength/velocity
        if self.velocity is None:
            self.__dict__.pop("velocity", None)
            if not isinstance(self.wavelength, Quantity):
                self.wavelength = self.wavelength * u.angstrom

        if self.wavelength is None:
            if not isinstance(self.velocity, Quantity):
                self.velocity = self.velocity * u.km_per_s
            self.__dict__.pop("wavelength", None)

        if not isinstance(self.flux, Quantity):
            self.flux = self.flux * u.Unit("1e-17 erg cm-2 s-1 Angstrom-1")

        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

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

    def vspec(self, wav_ref: float | Quantity, inplace: bool = True):
        from .methods.spectrum.radial_velocities import get_velocities

        struct = manage_inplace(self, inplace)

        if isinstance(wav_ref, Quantity):
            wav_ref = wav_ref.to(struct.wavelength.unit).value

        struct.velocity = get_velocities(struct.wavelength.value, wav_ref)
        struct.wav_ref = wav_ref * struct.wavelength.unit
        struct.wavelength = None

        return struct
