from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.Lightcurve import Lightcurve
    from ....structures.Powspec import Powspec

np.seterr(divide="ignore")


def gen_powspec(struct: object, lcs: list[Lightcurve], min: float, max: float, samples: int, multiband: bool = True) -> Powspec:
    from ....structures.Powspec import Powspec

    struct.kind = "powspec"

    if multiband:
        freq, power, fopt = do_ls(lcs, min, max, samples)
        band_str = ", ".join(list(set([lc.band for lc in lcs])))
        popt = (1 / fopt).to(u.day)

        return [
            Powspec(survey=struct.survey, band=band_str, frequency=freq, power=power, fopt=fopt, popt=popt, _target_key=lcs[0]._target_key)
        ]

    else:
        pspectra = []
        for lc in lcs:
            freq, power, fopt = do_ls(lc, min, max, samples)
            band_str = lc.band
            popt = (1 / fopt).to(u.day)

            pspectra.append(
                Powspec(survey=struct.survey, band=band_str, frequency=freq, power=power, fopt=fopt, popt=popt, _target_key=lc._target_key)
            )

        return pspectra
