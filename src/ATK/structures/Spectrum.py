from dataclasses import dataclass, field
from typing import Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .methods.spectrum.fitting import do_fitting
from .methods.spectrum.radial_velocities import get_rvs
from .structures_core import Container, manage_inplace
from .Target import Target


@dataclass(repr=False)
class Spectrum(Container):
    """
    Container for storing spectral data. This object stores both data and relevant metadata.

    .. rubric:: Valid Configurations

    A :class:`~ATK.Models.Spectrum` must be initialised with one of the following mutually exclusive parameters:

    - ``wavelength`` or ``velocity``

    Providing both or neither raises a ``ValueError``.
    """

    #: DOC_OVERRIDE
    survey: str | None = None
    #: DOC_OVERRIDE
    correction: str | None = None
    #: DOC_OVERRIDE
    search_pos: SkyCoord | None = None
    #: DOC_OVERRIDE
    separation: Quantity | None = None
    #: Exposure of spectrum.
    exposure: Quantity | None = None
    #: Reference wavelength.
    #:
    #: ``None`` unless spectrum has been converted to a velocity spectrum via :meth:`~ATK.Models.Spectrum.vspec`.
    wav_ref: Quantity | None = None

    #: Wavelength values.
    #: Mutually exclusive with ``velocity``.
    wavelength: Quantity | None = None
    #: Velocity values.
    #: Mutually exclusive with ``wavelength``.
    velocity: Quantity | None = None
    #: Flux values.
    flux: Quantity | None = None

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
    def _x_type(self):
        if self.wavelength is not None:
            return "wavelength"
        elif self.velocity is not None:
            return "velocity"
        else:
            raise ValueError("Spectrum container must hold one of 'wavelength' and 'velocity'.")

    @property
    def _x_arr(self):
        return getattr(self, self._x_type)

    def _set_x(self, val: numpy.ndarray):
        setattr(self, self._x_type, val)

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys = crop_nd(x=struct._x_arr, ys=ys, lower_lim=min, upper_lim=max)

        struct._set_x(x)
        struct.flux = ys[0]

        return struct

    crop.__doc__ = get_docstring("bin", x="``wavelength`` or ``velocity``", name="Spectrum")

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.binning import bin_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux]

        x, ys, _ = bin_nd(x=struct._x_arr, ys=ys, errs=[], bins=bins, size=size)

        struct._set_x(x)
        struct.flux = ys[0]

        return struct

    bin.__doc__ = get_docstring("bin", x="``wavelength`` or ``velocity``", name="Spectrum")

    def vspec(self, wav_ref: float | Quantity, inplace: bool = True):
        """
        Convert wavelength values to velocity relative to a reference wavelength.
        This replaces the stored ``wavelength`` axis with ``velocity``.

        Parameters
        ----------
        wav_ref : float | Quantity
            Reference (rest) wavelength used to compute velocities.

            If a :class:`~astropy.units.Unit` is not provided, ``wav_ref`` is assumed to be in the same unit as ``wavelength``.

        inplace : bool, optional
            If ``True``, modify the current :class:`~ATK.Models.Spectrum` in place - leaving the original unchanged. If ``False``, operate on and return a copy.
        """

        from .methods.spectrum.radial_velocities import get_velocities

        struct = manage_inplace(self, inplace)

        if isinstance(wav_ref, Quantity):
            wav_ref = wav_ref.to(struct.wavelength.unit).value

        struct.velocity = get_velocities(struct.wavelength.value, wav_ref)
        struct.wav_ref = wav_ref * struct.wavelength.unit
        struct.wavelength = None

        return struct

    @classmethod
    def from_dataframe(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_dataframe(target, data, **kwargs)

    from_dataframe.__func__.__doc__ = get_docstring("from_dataframe", obj="Spectrum", args=", ".join(f"``{p}``" for p in _required))

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        return super().from_table(target, data, **kwargs)

    from_table.__func__.__doc__ = get_docstring("from_table", obj="Spectrum", args=", ".join(f"``{p}``" for p in _required))
