from __future__ import annotations

from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from .timeseries_core import do_ls

if TYPE_CHECKING:
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


def compute_band_offsets(lcs, align_method="median"):
    ref_lc = lcs[0]

    if align_method == "mean":
        ref_val = np.average(ref_lc.brightness, weights=1 / np.power(ref_lc.brightness_err, 2))
    else:
        ref_val = np.median(ref_lc.brightness)

    band_offsets = {}

    for lc in lcs:
        if align_method == "mean":
            val = np.average(lc.brightness, weights=1 / np.power(lc.brightness_err, 2))
        else:
            val = np.median(lc.brightness)

        band_offsets[lc.band] = ref_val - val

    return band_offsets


def compute_phase_offset_midpoint(lcs, freq, align_method="median", nbins=50):
    all_phase = []
    all_flux = []

    for lc in lcs:
        phase = (lc.mjd * freq) % 1
        all_phase.append(phase)
        all_flux.append(lc.brightness)

    phase = np.concatenate(all_phase)
    flux = np.concatenate(all_flux)

    if align_method == "mean":
        midpoint = np.mean(flux)
    else:
        midpoint = np.median(flux)

    # binning
    bins = np.linspace(0, 1, nbins + 1)
    digitized = np.digitize(phase, bins)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.full(nbins, np.nan)

    for i in range(1, nbins + 1):
        vals = flux[digitized == i]
        if len(vals) > 0:
            bin_means[i - 1] = np.nanmean(vals)

    valid = ~np.isnan(bin_means)
    if not np.any(valid):
        return 0.0

    # choose bin closest to midpoint
    idx = np.where(valid)[0][np.argmin(np.abs(bin_means[valid] - midpoint))]
    return bin_centers[idx]


def compute_phase_offset_extrema(lcs, freq, mode="max"):
    all_phase = []
    all_flux = []

    brightness_type = list(set([lc.brightness_type for lc in lcs]))
    if len(brightness_type) > 1:
        raise ValueError("Invalid combination of lighcurve brightness types. Must be all 'flux' or all 'mag'.")
    brightness_type = brightness_type[0]

    for lc in lcs:
        phase = (lc.mjd * freq) % 1
        all_phase.append(phase)
        all_flux.append(lc.brightness)

    phase = np.concatenate(all_phase)
    flux = np.concatenate(all_flux)

    bins = np.linspace(0, 1, 50)
    digitized = np.digitize(phase, bins)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.array([np.nanmean(flux[digitized == i]) if np.any(digitized == i) else np.nan for i in range(1, len(bins))])

    valid = ~np.isnan(bin_means)
    phase = bin_centers[valid]
    flux = bin_means[valid]

    # sort
    idx = np.argsort(phase)
    phase = phase[idx]
    flux = flux[idx]

    # wrap phase for circular smoothing
    phase_ext = np.concatenate([phase - 1, phase, phase + 1])
    flux_ext = np.tile(flux, 3)

    # sort again (important after wrapping)
    idx = np.argsort(phase_ext)
    phase_ext = phase_ext[idx]
    flux_ext = flux_ext[idx]

    # smooth
    window = max(5, len(flux) // 20)
    kernel = np.ones(window) / window
    smooth_ext = np.convolve(flux_ext, kernel, mode="same")

    # restrict back to central region
    mask = (phase_ext >= 0) & (phase_ext < 1)

    phase = phase_ext[mask]
    smooth = smooth_ext[mask]

    # decide extrema direction
    if mode == "max":
        if brightness_type == "mag":
            i = np.argmin(smooth)
        else:
            i = np.argmax(smooth)
    elif mode == "min":
        if brightness_type == "mag":
            i = np.argmax(smooth)
        else:
            i = np.argmin(smooth)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return phase[i]


def compute_phase_offset(lcs, freq, align="median"):
    if align in ["mean", "median"]:
        return compute_phase_offset_midpoint(lcs, freq, align)
    elif align in ["max", "min"]:
        return compute_phase_offset_extrema(lcs, freq, align)
    else:
        raise ValueError(f"Unknown align mode: {align}")


def estimate_slope(phase, flux, window=0.1):
    # unwrap around 0 for continuity
    phase_ext = np.concatenate([phase - 1, phase, phase + 1])
    flux_ext = np.tile(flux, 3)

    mask = (phase_ext >= -window) & (phase_ext <= window)

    if np.sum(mask) < 5:
        return 0.0

    p = phase_ext[mask]
    f = flux_ext[mask]

    # simple linear fit
    coeffs = np.polyfit(p, f, 1)

    return coeffs[0]


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
    subtract: str | None = "median",
    repeat: int = 2,
    align: str = "median",
    freq: float | Quantity | None = None,
):
    from ....structures.Lightcurve import Lightcurve

    has_grid = (fmin is not None) and (fmax is not None) and (samples is not None)
    has_freq = freq is not None

    if not has_grid and not has_freq:
        raise ValueError("fold must be provided with (fmin, fmax, samples) or freq.")

    align_method = "mean" if subtract == "mean" else "median"

    fopts = {}
    phase_offsets = {}
    y_fits = {}

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

    if multiband:
        offset = compute_phase_offset(lcs, fopt.value, align)
        for lc in lcs:
            phase_offsets[lc.band] = offset
    else:
        for lc in lcs:
            fopt = fopts[lc.band]
            phase_offsets[lc.band] = compute_phase_offset([lc], fopt.value, align)

    if align == "mid":
        # choose reference band
        ref_band = lcs[0].band

        # compute reference slope
        ref_lc = next(lc for lc in lcs if lc.band == ref_band)
        fref = fopts[ref_band]

        phase_ref = (ref_lc.mjd * fref.value) % 1
        phase_ref = (phase_ref - phase_offsets[ref_band]) % 1

        slope_ref = estimate_slope(phase_ref, ref_lc.brightness)

        for lc in lcs:
            band = lc.band
            fopt = fopts[band]

            phase = (lc.mjd * fopt.value) % 1
            phase = (phase - phase_offsets[band]) % 1

            slope = estimate_slope(phase, lc.brightness)

            # flip if opposite orientation
            if slope * slope_ref < 0:
                phase_offsets[band] = (phase_offsets[band] + 0.5) % 1

    band_offsets = compute_band_offsets(lcs, align_method) if multiband else {}

    for band, y in y_fits.items():
        if y is None:
            continue

        if subtract == "mean":
            y_fits[band] -= np.mean(y)
        elif subtract == "median":
            y_fits[band] -= np.median(y)

    ctnrs = []
    for lc in lcs:
        fopt = fopts[lc.band]
        phase_offset = phase_offsets[lc.band]

        # phase
        phase = (lc.mjd * fopt.value) % 1
        phase = (phase - phase_offset) % 1

        # data
        brightness = lc.brightness

        if multiband:
            brightness = brightness + band_offsets.get(lc.band, 0)

        if subtract == "mean":
            weights = 1 / np.power(lc.brightness_err, 2)
            brightness = brightness - np.average(brightness, weights=weights)
        elif subtract == "median":
            brightness = brightness - np.median(brightness)

        phase, brightness, brightness_err = repeat_arrs(phase, brightness, lc.brightness_err, repeat=repeat)

        brightness_data = {f"{lc.brightness_type}": brightness, f"{lc.brightness_type}_err": brightness_err}

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
