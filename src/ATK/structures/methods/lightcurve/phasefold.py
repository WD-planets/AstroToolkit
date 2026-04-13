from __future__ import annotations

from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.Lightcurve import Lightcurve

np.seterr(divide="ignore")


def phase_dispersion_pdm(lc, freq, nbins=10, min_per_bin=5):
    """
    PDM theta statistic (lower = better period), Stellingwerf
    """

    # phase fold
    phase = (lc.mjd * freq) % 1
    flux = lc.brightness

    weights = 1 / (lc.brightness_err**2 + 1e-8)

    def weighted_var(x, w):
        mu = np.average(x, weights=w)
        return np.average((x - mu) ** 2, weights=w)

    global_var = weighted_var(flux, weights)

    # sort by phase
    idx = np.argsort(phase)
    phase = phase[idx]
    flux = flux[idx]

    # global variance (normalisation)
    global_var = np.var(flux)
    if global_var == 0:
        return np.inf

    # bin edges
    bins = np.linspace(0, 1, nbins + 1)
    digitized = np.digitize(phase, bins)

    # bin variance
    within_var = 0.0
    total_weight = 0

    for i in range(1, nbins + 1):
        mask = digitized == i

        vals = flux[mask]
        wvals = weights[mask]

        if len(vals) < min_per_bin:
            continue

        within_var += weighted_var(vals, wvals) * np.sum(wvals)
        total_weight += np.sum(wvals)

    if total_weight == 0:
        return np.inf

    within_var /= total_weight

    # theta statistic
    theta = within_var / global_var

    return theta


def optimise_freq(lcs, fopt, n_harmonics=5):
    if isinstance(fopt, Quantity):
        f0 = fopt.to_value(1 / u.day)
    else:
        f0 = fopt

    candidates = []

    for k in range(1, n_harmonics + 1):
        candidates.append(f0 * k)
        candidates.append(f0 / k)

    # remove duplicate / invalid candidate frequencies
    candidates = np.unique(np.array(candidates))
    candidates = candidates[candidates > 0]

    scores = []
    for f in candidates:
        scores.append(np.mean([phase_dispersion_pdm(lc, f) for lc in lcs]))

    best_freq = candidates[np.argmin(scores)]

    return best_freq * (1 / u.day)


def fold_lc(
    lcs: list[Lightcurve],
    fmin: float | Quantity | None = None,
    fmax: float | Quantity | None = None,
    samples: int | None = None,
    multiband: bool = True,
    optimise: bool = True,
    freq: float | Quantity | None = None,
):
    from ....structures.Lightcurve import Lightcurve

    has_grid = (fmin is not None) and (fmax is not None) and (samples is not None)
    has_freq = freq is not None

    if not has_grid and not has_freq:
        raise ValueError("fold must be provided with (fmin, fmax, samples) or freq.")

    fopts = {}

    if freq is not None:
        # user-specified frequency
        if isinstance(freq, Quantity):
            f_user = freq
        else:
            f_user = freq * (1 / u.day)

        for lc in lcs:
            fopt = f_user
            fopts[lc.band] = f_user

    else:
        if multiband:
            _, _, fopt, ls = do_ls(lcs, fmin, fmax, samples, return_model=True)
            if optimise:
                fopt = optimise_freq(lcs, fopt)

            for lc in lcs:
                fopts[lc.band] = fopt

        else:
            ls_models = {}

            for lc in lcs:
                _, _, fopt, ls = do_ls(lc, fmin, fmax, samples, return_model=True)
                if optimise:
                    fopt = optimise_freq([lc], fopt)

                fopts[lc.band] = fopt
                ls_models[lc.band] = ls

    ctnrs = []
    for lc in lcs:
        fopt = fopts[lc.band]

        # phase
        phase = (lc.mjd * fopt.value) % 1

        # data
        brightness = lc.brightness

        brightness_data = {f"{lc.brightness_type}": brightness, f"{lc.brightness_type}_err": lc.brightness_err}

        f_lc = Lightcurve(
            survey=lc.survey,
            band=lc.band,
            obj_id=lc.obj_id,
            multiband=multiband,
            phase=phase,
            fopt=fopts[lc.band],
            popt=1 / fopts[lc.band],
            _target_key=lc._target_key,
            **brightness_data,
        )

        ctnrs.append(f_lc)

    return ctnrs
