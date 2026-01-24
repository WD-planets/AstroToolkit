from __future__ import annotations

from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
import pandas as pd

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.definitions import Lightcurve


def fold_lc(struct: object, lcs: list[Lightcurve], min: float, max: float, samples: int):
    from ....structures.definitions import FoldedLightcurve

    freq, power, fopt, ls = do_ls(lcs, min, max, samples, return_model=True)

    ctnrs = []
    for lc in lcs:
        t_fit = np.linspace(0, 1 / fopt.value, 1000) * u.day
        phase = lc.mjd % (1 / fopt.to(1 / u.day).value)
        y_fit = ls.model(t=t_fit, frequency=fopt)

        # data = pd.DataFrame({"line_x": t_fit, "line_y": y_fit})

        ms_brightness = lc.brightness - np.median(lc.brightness)
        brightness_err = lc.brightness_err

        brightness_data = {f"ms_{lc.brightness_type}": ms_brightness, f"ms_{lc.brightness_type}_err": brightness_err}

        f_lc = FoldedLightcurve(survey=lc.survey, band=lc.band, obj_id=lc.obj_id, phase=phase, fopt=fopt, popt=1 / fopt, **brightness_data)

        ctnrs.append(f_lc)

    return ctnrs
