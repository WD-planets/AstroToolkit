from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.timeseries import LombScargle, LombScargleMultiband

# get around circular import for type hinting
if TYPE_CHECKING:
    from ....structures.Lightcurve import Lightcurve

warnings.filterwarnings("ignore", category=RuntimeWarning)


def format_lc_data(lcs: list[Lightcurve] | Lightcurve):
    if isinstance(lcs, list):
        mjd, brightness, brightness_err, band = [], [], [], []
        for lc in lcs:
            mjd.append(lc.mjd)
            brightness.append(lc.brightness)
            brightness_err.append(lc.brightness_err)
            band.append(np.full(len(lc.mjd), lc.band))
    else:
        mjd, brightness, brightness_err, band = (lcs.mjd, lcs.brightness, lcs.brightness_err, np.full(len(lcs.mjd), lcs.band))

    mjd = np.concatenate(mjd, axis=None) * u.day
    brightness = np.concatenate(brightness, axis=None)
    brightness_err = np.concatenate(brightness_err, axis=None)
    band = np.concatenate(band, axis=None)

    return mjd, brightness, brightness_err, band


def do_ls(lcs: list[Lightcurve] | Lightcurve, min: float, max: float, samples: int, return_model: bool = False):
    mjd, brightness, brightness_err, band = format_lc_data(lcs)

    freqs = np.linspace(min, max, samples) * (1 / u.day)
    if isinstance(lcs, list):
        ls = LombScargleMultiband(t=mjd, y=brightness, dy=brightness_err, bands=band, nterms_base=1, nterms_band=0)
    else:
        ls = LombScargle(mjd, brightness, brightness_err)

    power = ls.power(freqs).value

    if np.isnan(power).all():
        best_freq = 0.0 * (1 / u.day)
    else:
        best_freq = freqs[np.nanargmax(power)]

    if return_model:
        return freqs, power, best_freq, ls

    return freqs, power, best_freq
