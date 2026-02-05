from __future__ import annotations

from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.DataSet import DataSet
    from ....structures.Lightcurve import Lightcurve

np.seterr(divide="ignore")


def repeat_arrs(x: np.ndarray, *ys: np.ndarray, repeat: int = 2):
    xs = []
    ys_out = [[] for _ in ys]

    for k in range(repeat):
        xs.append(x + k)
        for i, y in enumerate(ys):
            ys_out[i].append(y)

    x_rep = np.concatenate(xs)
    ys_rep = tuple(np.concatenate(ylist) for ylist in ys_out)

    return (x_rep, *ys_rep)


def fold_lc(
    lcs: list[Lightcurve],
    min: float,
    max: float,
    samples: int,
    multiband: bool = True,
    subtract: str | None = None,
    fit: bool = False,
    repeat: int = 2,
):
    from ....structures.Lightcurve import Lightcurve

    align_method = "mean" if subtract == "mean" else "median"

    phase_fit = np.linspace(0, 1, 100)
    if multiband:
        y_fits = {}
        freq, power, fopt, ls = do_ls(lcs, min, max, samples, return_model=True)
        t_fit = (phase_fit / fopt.to(1 / u.day).value) * u.day
        y_fit = ls.model(t_fit, fopt)
        y_fits = {str(band): y_fit[i] for i, band in enumerate(set(ls.bands))}
    else:
        y_fits = {}
        for lc in lcs:
            freq, power, fopt, ls = do_ls(lc, min, max, samples, return_model=True)
            t_fit = (phase_fit / fopt.to(1 / u.day).value) * u.day
            try:
                y_fit = ls.model(t_fit, fopt)
            except np.linalg.LinAlgError:
                y_fit = None
            y_fits[lc.band] = y_fit

    # offsets only needed if using MultiBand
    if multiband:
        # Compute per-band photometry weighted mean
        if align_method == "mean":
            band_centroids = {}
            for lc in lcs:
                weights = 1 / np.pow(lc.brightness_err, 2)
                band_centroids[lc.band] = np.average(lc.brightness, weights=weights)
        elif align_method == "median":
            band_centroids = {lc.band: np.median(lc.brightness) for lc in lcs}

        # Compute per-band sinusoid mean from LS model
        if align_method == "mean":
            sinusoid_centroids = {band: np.mean(y_fits[band]) for band in y_fits}
        elif align_method == "median":
            sinusoid_centroids = {band: np.median(y_fits[band]) for band in y_fits}

        # Compute the offset needed to match photometry
        band_offsets = {band: band_centroids[band] - sinusoid_centroids[band] for band in y_fits}

        for band in y_fits:
            if y_fits[band] is not None:
                y_fits[band] += band_offsets[band]

    for band in y_fits:
        if subtract == "mean" and y_fits[band] is not None:
            y_fits[band] = y_fits[band] - np.mean(y_fits[band])
        elif subtract == "median" and y_fits[band] is not None:
            y_fits[band] = y_fits[band] - np.median(y_fits[band])

    phase_offset = phase_fit[np.nanargmax(y_fits[list(y_fits.keys())[0]])]

    ctnrs = []
    for lc in lcs:
        phase = (lc.mjd * fopt.to(1 / u.day).value) % 1

        phase = (phase - phase_offset) % 1
        fit_x = (np.linspace(0, 1, 100) - phase_offset) % 1
        order = np.argsort(fit_x)

        if y_fits[lc.band] is not None:
            fit_x = fit_x[order]
            fit_y = y_fits[lc.band][order]
        else:
            fit_x = None
            fit_y = None

        if subtract == "mean":
            weights = 1 / np.pow(lc.brightness_err, 2)
            brightness = lc.brightness - np.average(lc.brightness, weights=weights)
        elif subtract == "median":
            brightness = lc.brightness - np.median(lc.brightness)
        else:
            brightness = lc.brightness

        if y_fits[lc.band] is not None:
            fit_x, fit_y = repeat_arrs(fit_x, fit_y, repeat=2)
        phase, brightness, brightness_err = repeat_arrs(phase, brightness, lc.brightness_err, repeat=repeat)

        brightness_data = {f"{lc.brightness_type}": brightness, f"{lc.brightness_type}_err": brightness_err}

        f_lc = Lightcurve(
            survey=lc.survey,
            band=lc.band,
            obj_id=lc.obj_id,
            phase=phase,
            fit_x=fit_x,
            fit_y=fit_y,
            fopt=fopt,
            popt=1 / fopt,
            _target_key=lc._target_key,
            **brightness_data,
        )

        if not fit:
            f_lc.fit_x = None
            f_lc.fit_y = None

        ctnrs.append(f_lc)

    return ctnrs
