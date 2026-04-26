from dataclasses import dataclass, field
from typing import ClassVar, Self

import astropy.units as u
import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .Base import (Container, DataFrameIOMixin, FITSIOMixin, TableIOMixin,
                   manage_inplace)
from .methods.lightcurve.phasefold import fold_lc
from .methods.lightcurve.powspec import gen_powspec
from .Powspec import Powspec
from .Target import Target

PLOT_PARAMS = {
    "bands": [
        "list of str, optional",
        """
        Bands to process and include in figure. Only used 
        
        If not provided, all bands are plotted. See :doc:`here </auto_tutorials/lightcurves/lightcurve_plotting>` for a list of supported bands in default surveys.
        """,
    ],
    "colours": [
        "list of str, optional",
        """
    Colour to give to each band.

    If not provided, a default set of colours are used. The following colours are supported: ``"green"``, ``"red"``, ``"blue"``, ``"orange"``, ``"purple"``, ``"black"``. 

    If ``bands`` is ``None``, colours are automatically assigned among all available bands.
    """,
    ],
    "cmap": [
        "{'mean', 'flat'}, optional",
        """
        Sets the colour map. Only used in unfolded :class:`~ATK.Models.Lightcurve`\ s, i.e. when ``phase`` is ``None``.

        ``'mean'``: colour scales with distance from the centre.

        ``'flat'``: no colour scaling.

        Default is ``'mean'``.
    """,
    ],
    "time_format": [
        "{'reduced', 'original'}, optional",
        """
        Sets the format of the time axis. Only used in unfolded :class:`~ATK.Models.Lightcurve`\ s, i.e. when ``phase`` is ``None``.

        ``'reduced'``: minimum MJD across all available photometry is subtracted. Time axis starts at date of first observation.

        ``'original'``: time axis is MJD.

        Default is ``'reduced'``.
        """,
    ],
    "subtract": [
        "{'mean', 'mean', None}, optional",
        """
        Sets metric used in subtracting (and hence aligning) photometry. E.g. if ``'median``, y-axis represents change in brightness relative to median. Only used in folded :class:`~ATK.Models.Lightcurve`\ s, i.e. when ``phase`` is not ``None``.

        If ``None``, no photometry subtraction/alignment is performed.
         
        Default is ``'median'``.
        """,
    ],
    "align": [
        "{'max', 'min', 'median', 'mean', None}, optional",
        """
        Aligns photometry to a chosen feature.

        ``'max'``: photometry is aligned to begin at a local maximum in modulation.

        ``'min'``: photometry is aligned to begin at a local minimum in modulation.

        ``'median'``: photometry is aligned to begin at the median.

        ``'mean'``: photometry is aligned to begin at the mean.
        
        If ``None``, no phase-alignment is performed.

        Default is ``'max'``.
        """,
    ],
    "repeat": [
        "int, optional",
        """
        Sets the number of modulations in brightness to display.

        Default is ``2``.
        """,
    ],
}


@dataclass(repr=False)
class Lightcurve(Container, DataFrameIOMixin, TableIOMixin, FITSIOMixin):
    """
    Container for storing time-series photometry. This object stores both data and relevant metadata.

    |

    .. rubric:: Valid Configurations

    A :class:`~ATK.Models.Lightcurve` must be initialised with one of the following mutually exclusive forms:

    **Photometry:**

    - ``flux`` and ``flux_err``
    - ``mag`` and ``mag_err``

    **Time axis:**

    - ``mjd``
    - ``phase``

    Providing both or neither in either group raises a ``ValueError``.

    |

    """

    #: DOC_OVERRIDE
    survey: str | None = None
    #: DOC_OVERRIDE
    band: str | None = None
    #: DOC_OVERRIDE
    correction: str | None = None
    #: DOC_OVERRIDE
    search_pos: SkyCoord | None = None
    #: DOC_OVERRIDE
    separation: Quantity | None = None
    #: Survey-specific object ID.
    #:
    #: Set when performing light curve queries with ``split = True``.
    obj_id: str | None = None
    _multiband: bool | None = None

    _data_methods: tuple = ("crop", "bin", "clip")
    _group_data_methods: dict = field(default_factory=lambda: {"fold": fold_lc, "pspec": gen_powspec})
    _group_data_methods_doc: ClassVar[dict] = {"fold": fold_lc, "pspec": gen_powspec}

    #: Modified Julian Day values.
    #: Mutually exclusive with ``phase``.
    mjd: Quantity | None = None
    #: Flux values.
    #: Mutually exclusive with ``mag``.
    flux: Quantity | None = None
    #: Flux error values.
    #: Mutually exclusive with ``mag_err``.
    flux_err: Quantity | None = None
    #: Magnitude values.
    #: Mutually exclusive with ``flux``.
    mag: Quantity | None = None
    #: Magnitude error values.
    #: Mutually exclusive with ``flux_err``.
    mag_err: Quantity | None = None
    #: Right ascension values.
    #:
    #: ``None`` unless light curve has been folded with :meth:`~ATK.Models.Lightcurve.fold`.
    ra: Quantity | None = None
    #: Declination values.
    #:
    #: ``None`` unless light curve has been folded with :meth:`~ATK.Models.Lightcurve.fold`.
    dec: Quantity | None = None

    #: Phase values.
    #: Mutually exclusive with ``mjd``.
    phase: Quantity | None = None
    #: Fold frequency.
    #:
    #: Only relevant in folded light curves (i.e. when ``phase`` is not ``None``).
    fopt: Quantity | None = None
    #: Fold period.
    #:
    #: Only relevant in folded light curves (i.e. when ``phase`` is not ``None``).
    popt: Quantity | None = None

    _required = ["survey", "band"]

    _units = {
        "mjd": u.day,
        "mag": u.mag,
        "mag_err": u.mag,
        "flux": u.count / u.s,
        "flux_err": u.count / u.s,
        "ra": u.deg,
        "dec": u.deg,
        "phase": u.one,
    }

    _plot_params = PLOT_PARAMS

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

    def crop(self, cmin: float | Quantity | None = None, cmax: float | Quantity | None = None, inplace: bool = True) -> Self:
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct._brightness, struct._brightness_err]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                ys.append(val)

        x, ys = crop_nd(x=struct._time, ys=ys, lower_lim=cmin, upper_lim=cmax)

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

    crop.__doc__ = get_docstring("crop", x="``mjd`` or ``phase``", name="Lightcurve")

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True) -> Self:
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

    bin.__doc__ = get_docstring("bin", x="``mjd`` or ``phase``", name="Lightcurve")

    def clip(self, sigma: float, sigma_lower: float | None = None, sigma_upper: float | None = None, inplace: bool = False) -> Self:
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

    clip.__doc__ = get_docstring("clip", y="``mag`` or ``flux``", name="Lightcurve")

    def pspec(self, fmin: float | Quantity, fmax: float | Quantity, samples: int) -> Powspec:
        """
        Generates a :class:`~ATK.Models.Powspec` (power spectrum) using :class:`astropy.timeseries.LombScargle` over a range of trial frequencies.

        Parameters
        ----------
        fmin : float | Quantity
            Minimum frequency.

            If a :class:`~astropy.units.Unit` is not provided, ``fmin`` is assumed to be in :math:`\mathrm{days}^{-1}`.

        fmax : float | Quantity
            Maximum frequency.

            If a :class:`~astropy.units.Unit` is not provided, ``fmax`` is assumed to be in :math:`\mathrm{days}^{-1}`.

        samples: int
            Number of samples in frequency range defined by ``fmin`` and ``fmax``.

        Returns
        -------
        :class:`~ATK.Models.Powspec`
        """

        struct = gen_powspec([self], fmin, fmax, samples)[0]

        return struct

    def fold(
        self,
        fmin: float | Quantity,
        fmax: float | Quantity,
        samples: int,
        optimise: bool = True,
        freq: float | Quantity | None = None,
        inplace: bool = True,
    ) -> Self:
        """
        Phase-folds the :class:`~ATK.Models.Lightcurve` curve onto the best frequency over a range of trial frequencies, or a fixed frequency.

        Parameters
        ----------
        fmin : float or :class:`~astropy.units.Quantity`, optional
            Minimum frequency.

            If a :class:`~astropy.units.Quantity` is not provided, ``fmin`` is assumed
            to be in :math:`\mathrm{days}^{-1}`.

        fmax : float or :class:`~astropy.units.Quantity`, optional
            Maximum frequency.

            If a :class:`~astropy.units.Quantity` is not provided, ``fmax`` is assumed
            to be in :math:`\mathrm{days}^{-1}`.

        samples : int, optional
            Number of samples in the frequency range defined by ``fmin`` and ``fmax``.

        optimise : bool, optional
            If True, a phase-dispersion metric is calculated for a set of harmonics either side of the peak frequency, and the best is chosen. This can help to preserve real periodic structure, especially in the case of asymmetric modulation.

        freq : float or :class:`~astropy.units.Quantity`, optional
            Specific frequency on which to fold the :class:`~ATK.Models.Lightcurve`. If provided, this overrides ``fmin``, ``fmax``, and ``samples``.

            If a :class:`~astropy.units.Quantity` is not provided, ``freq`` is assumed
            to be in :math:`\mathrm{days}^{-1}`.

        inplace : bool, optional
            If ``True``, modify the current :class:`~ATK.Models.Lightcurve` in place - leaving the original unchanged. If ``False``, operate on and return a copy.

        Returns
        -------
        ``Self``
            The folded :class:`~ATK.Models.Lightcurve`, with ``fopt`` and ``popt`` set to the folded frequency and folded period, respectively. Returns ``self`` if ``inplace=True``, otherwise returns a new instance.
        """

        struct = manage_inplace(self, inplace)
        struct = fold_lc([struct], fmin=fmin, fmax=fmax, samples=samples, multiband=False, optimise=optimise, freq=freq)[0]

        # mutate self to match folded lightcurve, as fold_lc creates a new object
        if inplace:
            self.__dict__.clear()
            self.__dict__.update(struct.__dict__)

            return self

        return struct
