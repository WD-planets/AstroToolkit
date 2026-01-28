from dataclasses import dataclass, field

import numpy
from astropy.coordinates import SkyCoord
from astropy.units import Quantity, u

from .structures_core import (GROUP_METHODS, Container, QuantityArray,
                              manage_inplace)


@dataclass(repr=False)
class Lightcurve(Container):
    # --- metadata ---
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    separation: Quantity | None = None
    band: str | None = None
    obj_id: str | None = None

    _required: list = field(default_factory=list)

    # --- data ---
    mjd: numpy.ndarray | None = None
    flux: numpy.ndarray | QuantityArray | None = None
    flux_err: numpy.ndarray | QuantityArray | None = None
    mag: numpy.ndarray | None = None
    mag_err: numpy.ndarray | None = None
    ra: numpy.ndarray | None = None
    dec: numpy.ndarray | None = None

    # --- folded data ---
    phase: numpy.ndarray | None = None
    fit_x: numpy.ndarray | None = None
    fit_y: numpy.ndarray | None = None
    fopt: Quantity | None = None
    popt: Quantity | None = None

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def __post_init__(self):
        # check for a valid input combination
        if (self.flux is None) == (self.mag is None):
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")

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

    @property
    def brightness(self):
        return getattr(self, self.brightness_type)

    @property
    def brightness_err(self):
        return getattr(self, f"{self.brightness_type}_err")

    def set_brightness(self, val: numpy.ndarray):
        setattr(self, self.brightness_type, val)

    def set_brightness_err(self, val: numpy.ndarray):
        setattr(self, f"{self.brightness_type}_err", val)

    @property
    def brightness_type(self):
        if self.mag is not None:
            return "mag"
        elif self.flux is not None:
            return "flux"
        else:
            # shouldn't happen due to __post_init__
            raise ValueError("Lightcurve container must hold one of 'mag' and 'flux'.")

    @property
    def time_type(self):
        if self.mjd is not None:
            return "mjd"
        elif self.phase is not None:
            return "phase"
        else:
            raise ValueError("Lightcurve container must hold one of 'mjd' and 'phase'.")

    @property
    def time(self):
        return getattr(self, self.time_type)

    def set_time(self, val: numpy.ndarray):
        setattr(self, self.time_type, val)

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.brightness, struct.brightness_err]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                ys.append(val)

        x, ys = crop_nd(x=struct.time, ys=ys, lower_lim=min, upper_lim=max)

        if len(ys) > 2:
            brightness, brightness_err, ra, dec = ys
        else:
            brightness, brightness_err = ys

        struct.set_time(x)
        struct.set_brightness(brightness)
        struct.set_brightness_err(brightness_err)

        if getattr(struct, "fit_x", None) is not None and getattr(struct, "fit_x", None) is not None:
            ys = [struct.fit_y]
            x, ys = crop_nd(struct.fit_x, ys=ys, lower_lim=min, upper_lim=max)

            struct.fit_x = x
            struct.fit_y = ys[0]

        return struct

    def bin(self, bins: int | None = None, size: Quantity | float | None = None, inplace=True):
        from .methods.binning import bin_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.brightness]
        for attr in ["ra", "dec"]:
            val = getattr(struct, attr)
            if val is not None:
                ys.append(val)

        if size is not None and not isinstance(size, Quantity):
            size = size * u.day

        x, ys, errs = bin_nd(x=struct.time, ys=ys, errs=[struct.brightness_err], bins=bins, size=size)

        if len(ys) > 1:
            brightness, ra, dec = ys
            struct.ra = ra
            struct.dec = dec
        else:
            brightness = ys[0]
        brightness_err = errs[0]

        struct.set_time(x)
        struct.set_brightness(brightness)
        struct.set_brightness_err(brightness_err)

        return struct

    def fold(ctnrs: list[object], min: float, max: float, samples: int):
        return GROUP_METHODS["fold"](ctnrs, min=min, max=max, samples=samples)

    def pspec(self, samples: int): ...
