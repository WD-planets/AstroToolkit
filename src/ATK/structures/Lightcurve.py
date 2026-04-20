from dataclasses import dataclass, field

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from ..utilities.docstrings import get_docstring
from .methods.lightcurve.phasefold import fold_lc
from .methods.lightcurve.powspec import gen_powspec
from .structures_core import Container, manage_inplace


@dataclass(repr=False)
class Lightcurve(Container):
    """
    Container for storing time-series photometry. This object stores both data and relevant metadata. test

    .. rubric:: Valid Combinations

    A :class:`~ATK.Models.Lightcurve` must be initialized with one of the following mutually exclusive forms:

    **Photometry:**

    - ``flux`` and ``flux_err``
    - ``mag`` and ``mag_err``

    **Time axis:**

    - ``mjd``
    - ``phase``

    Providing both or neither in either group raises a ``ValueError``.

    |

    """

    #: Survey from which light curve originates.
    survey: str | None = None
    #: Photometric band of stored data.
    band: str | None = None
    #: Achieved degree of proper motion correction.
    #:
    #: - ``'full'`` = complete 3-dimensional projection on the sky.
    #: - ``'partial'`` = 2-dimensional plane projection.
    #: - ``'none'`` = no correction.
    correction: str | None = None
    #: Position of search at time of execution (i.e. post-correction).
    search_pos: SkyCoord | None = None
    #: Separation between position of the search and the returned data.
    separation: Quantity | None = None
    #: Survey-specific object ID.
    #:
    #: Set when performing light curve queries with ``split = True``.
    obj_id: str | None = None
    _multiband: bool | None = None

    _data_methods: tuple = ("crop", "bin", "clip")
    _group_data_methods: dict = field(default_factory=lambda: {"fold": fold_lc, "pspec": gen_powspec})

    #: Modified Julian Day values.
    #: Mutually exclusive with ``phase``.
    mjd: numpy.ndarray | None = None
    #: Flux values.
    #: Mutually exclusive with ``mag``.
    flux: numpy.ndarray | Quantity | None = None
    #: Flux error values.
    #: Mutually exclusive with ``mag_err``.
    flux_err: numpy.ndarray | Quantity | None = None
    #: Magnitude values.
    #: Mutually exclusive with ``flux``.
    mag: numpy.ndarray | None = None
    #: Magnitude error values.
    #: Mutually exclusive with ``flux_err``.
    mag_err: numpy.ndarray | None = None
    #: Right ascension values.
    ra: numpy.ndarray | None = None
    #: Declination values.
    dec: numpy.ndarray | None = None

    #: Phase values.
    #: Mutually exclusive with ``mjd``.
    phase: numpy.ndarray | None = None
    #: Fold frequency.
    #:
    #: Only relevant in folded light curves (i.e. when ``phase`` is not ``None``).
    fopt: Quantity | None = None
    #: Fold period.
    #:
    #: Only relevant in folded light curves (i.e. when ``phase`` is not ``None``).
    popt: Quantity | None = None

    _required = ["survey", "band"]

    _units = {"mjd": u.day, "mag": u.mag, "mag_err": u.mag, "flux": u.count / u.s, "flux_err": u.count / u.s, "ra": u.deg, "dec": u.deg}

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def __post_init__(self):
        # check for a valid input combination
        if (self.flux is None) == (self.mag is None):
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")

        if self.flux is not None:
            if self.flux_err is None:
                raise ValueError("Lightcurve container is missing flux_err")
        if self.mag is not None:
            if self.mag_err is None:
                raise ValueError("Lightcurve container is missing mag_err")

        # ensure valid combination of flux/flux_err/mag/mag_err
        if self.flux is None:
            for f in ("flux", "flux_err"):
                self.__dict__.pop(f, None)
        if self.mag is None:
            for f in ("mag", "mag_err"):
                self.__dict__.pop(f, None)

        if (self.mjd is None) == (self.phase is None):
            raise ValueError("Lightcurve container must hold one of 'mjd' and 'phase'.")

        if self.phase is None:
            self.__dict__.pop("phase", None)
        if self.mjd is None:
            self.__dict__.pop("mjd", None)

        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    @property
    def _brightness(self):
        return getattr(self, self._brightness_type)

    @property
    def _brightness_err(self):
        return getattr(self, f"{self._brightness_type}_err")

    def _set_brightness(self, val: numpy.ndarray):
        setattr(self, self._brightness_type, val)

    def _set_brightness_err(self, val: numpy.ndarray):
        setattr(self, f"{self._brightness_type}_err", val)

    @property
    def _brightness_type(self):
        if self.mag is not None:
            return "mag"
        elif self.flux is not None:
            return "flux"
        else:
            # shouldn't happen due to __post_init__
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")

    @property
    def _time_type(self):
        if self.mjd is not None:
            return "mjd"
        elif self.phase is not None:
            return "phase"
        else:
            raise ValueError("Lightcurve container must hold one of 'mjd' and 'phase'.")

    @property
    def _time(self):
        return getattr(self, self._time_type)

    def _set_time(self, val: numpy.ndarray):
        setattr(self, self._time_type, val)

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct._brightness, struct._brightness_err]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                ys.append(val)

        x, ys = crop_nd(x=struct._time, ys=ys, lower_lim=min, upper_lim=max)

        if len(ys) > 2:
            brightness, brightness_err, ra, dec = ys
            struct.ra = ra
            struct.dec = dec
        else:
            brightness, brightness_err = ys

        struct._set_time(x)
        struct._set_brightness(brightness)
        struct._set_brightness_err(brightness_err)

        return struct

    crop.__doc__ = get_docstring("crop", x="mjd", name="Lightcurve")

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.binning import bin_nd

        struct = manage_inplace(self, inplace)

        ys = [struct._brightness]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                ys.append(val)

        x, ys, errs = bin_nd(x=struct._time, ys=ys, errs=[struct._brightness_err], bins=bins, size=size)

        if len(ys) > 1:
            brightness, ra, dec = ys
            struct.ra = ra
            struct.dec = dec
        else:
            brightness = ys[0]
        brightness_err = errs[0]

        struct._set_time(x)
        struct._set_brightness(brightness)
        struct._set_brightness_err(brightness_err)

        return struct

    bin.__doc__ = get_docstring("bin", x="mjd", name="Lightcurve")

    def clip(self, sigma: float, sigma_lower: float | None = None, sigma_upper: float | None = None, inplace: bool = False):
        from .methods.sigma_clip import do_sigma_clipping

        struct = manage_inplace(self, inplace)

        arrs = [struct.mjd, struct._brightness_err]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                arrs.append(val)

        y, arrs = do_sigma_clipping(y=struct._brightness, arrs=arrs, sigma=sigma, sigma_lower=sigma_lower, sigma_upper=sigma_upper)

        if len(arrs) > 2:
            mjd, brightness_err, ra, dec = arrs
            struct.ra = ra
            struct.dec = dec
        else:
            mjd, brightness_err = arrs

        struct._set_time(mjd)
        struct._set_brightness(y)
        struct._set_brightness_err(brightness_err)

        return struct

    clip.__doc__ = get_docstring("clip", name="Lightcurve")

    def pspec(self, fmin: float, fmax: float, samples: int):
        """
        Generates a :class:`~ATK.Models.Powspec` using :class:`astropy.timeseries.LombScargle`.

        Parameters
        ----------
        fmin : float
        """

        struct = gen_powspec([self], fmin, fmax, samples)[0]

        return struct

    def fold(
        self, fmin: float, fmax: float, samples: int, optimise: bool = True, freq: float | Quantity | None = None, inplace: bool = True
    ):
        struct = manage_inplace(self, inplace)
        struct = fold_lc([struct], fmin=fmin, fmax=fmax, samples=samples, multiband=False, optimise=optimise, freq=freq)[0]

        # mutate self to match folded lightcurve, as fold_lc creates a new object
        if inplace:
            self.__dict__.clear()
            self.__dict__.update(struct.__dict__)

            return self

        return struct
