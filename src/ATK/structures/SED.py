from dataclasses import dataclass, field

import numpy
from astropy.units import Quantity

from .structures_core import Container, manage_inplace


@dataclass(repr=False)
class SED(Container):
    survey: numpy.ndarray | None = None
    correction: numpy.ndarray | None = None
    band: numpy.ndarray | None = None
    id: numpy.ndarray | None = None
    separation: numpy.ndarray | Quantity | None = None
    wavelength: numpy.ndarray | Quantity | None = None
    flux: numpy.ndarray | Quantity | None = None
    flux_err: numpy.ndarray | Quantity | None = None

    _data_methods: tuple = ("crop", "bin")

    def __repr__(self):
        return "<Spectral Energy Distribution>"

    def crop(self, min: float | None = None, max: float | None = None, inplace=True):
        from .methods.cropping import crop_nd

        struct = manage_inplace(self, inplace)

        ys = [struct.flux, struct.flux_err, struct.survey, struct.correction, struct.band, struct.separation]

        x, ys = crop_nd(x=struct.wavelength, ys=ys, lower_lim=min, upper_lim=max)

        flux, flux_err, survey, correction, band, separation = ys

        struct.wavelength = x
        struct.flux = flux
        struct.flux_err = flux_err
        struct.survey = survey
        struct.correction = correction
        struct.band = band
        struct.separation = separation

        return struct
