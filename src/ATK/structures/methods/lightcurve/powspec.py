from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u

from ....utilities.units import _align_to_unit, _strip_unit
from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.Lightcurve import Lightcurve
    from ....structures.Powspec import Powspec

np.seterr(divide="ignore")


def gen_powspec(lcs: list[Lightcurve], fmin: float, fmax: float, samples: int, multiband: bool = True) -> Powspec:
    from ....structures.Powspec import Powspec

    _, mjd_unit = _strip_unit(lcs[0].mjd)
    freq_unit = (1 / mjd_unit) if mjd_unit is not None else None
    fmin = _align_to_unit(fmin, freq_unit, "fmin", "lcs.mjd")
    fmax = _align_to_unit(fmax, freq_unit, "fmax", "lcs.mjd")

    if multiband:
        freq, power, fopt = do_ls(lcs, fmin, fmax, samples)
        band_str = ", ".join(list(set([lc.band for lc in lcs])))
        popt = (1 / fopt).to(u.day)

        return [
            Powspec(survey=lcs[0].survey, band=band_str, frequency=freq, power=power, fopt=fopt, popt=popt, _target_key=lcs[0]._target_key)
        ]

    else:
        pspectra = []
        for lc in lcs:
            freq, power, fopt = do_ls(lc, fmin, fmax, samples)
            band_str = lc.band
            popt = (1 / fopt).to(u.day)

            pspectra.append(
                Powspec(
                    survey=lcs[0].survey,
                    obj_id=lc.obj_id,
                    band=band_str,
                    _multiband=multiband,
                    frequency=freq,
                    power=power,
                    fopt=fopt,
                    popt=popt,
                    _target_key=lc._target_key,
                )
            )

        return pspectra
