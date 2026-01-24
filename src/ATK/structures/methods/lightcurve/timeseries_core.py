from __future__ import annotations

from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.timeseries import LombScargleMultiband

# get around circular import for type hinting
if TYPE_CHECKING:
    from ....structures.definitions import Lightcurve


def format_lc_data(lcs: list[Lightcurve]):
    mjd, brightness, brightness_err, band = [], [], [], []

    for lc in lcs:
        mjd.append(lc.mjd)
        brightness.append(lc.brightness)
        brightness_err.append(lc.brightness_err)
        band.append(np.full(len(lc.mjd), lc.band))

    mjd = np.concatenate(mjd, axis=None) * u.day
    brightness = np.concatenate(brightness, axis=None)
    brightness_err = np.concatenate(brightness_err, axis=None)
    band = np.concatenate(band, axis=None)

    return mjd, brightness, brightness_err, band


def do_ls(lcs: list[Lightcurve], min: float, max: float, samples: int, return_model: bool):
    mjd, brightness, brightness_err, band = format_lc_data(lcs)

    freqs = np.linspace(min, max, samples) * (1 / u.day)
    ls = LombScargleMultiband(mjd, brightness, band, brightness_err, nterms_base=1, nterms_band=0)
    power = ls.power(freqs)

    if np.isnan(power).all():
        return None

    best_freq = freqs[np.argmax(power)]

    if return_model:
        return freqs, power, best_freq, ls

    return freqs, power, best_freq
