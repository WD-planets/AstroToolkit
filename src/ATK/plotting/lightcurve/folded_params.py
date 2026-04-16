import copy

import numpy as np

from ...structures.Lightcurve import Lightcurve


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

    if align_method is None:
        band_offsets = {}
        for lc in lcs:
            band_offsets[lc.band] = 0.0
        return band_offsets

    if align_method not in ["mean", "median"]:
        raise ValueError(f"Unexpected align method '{align_method}'. Accepted methods: 'mean', 'median'.")

    if align_method == "mean":
        ref_val = np.average(ref_lc.brightness, weights=1 / np.power(ref_lc.brightness_err, 2))
    elif align_method == "median":
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
        phase = lc.phase
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
        phase = lc.phase
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


def dispatch_phase_offsets(lcs, freq, align="median"):
    if align in ["mean", "median"]:
        return compute_phase_offset_midpoint(lcs, freq, align)
    elif align in ["max", "min"]:
        return compute_phase_offset_extrema(lcs, freq, align)
    else:
        raise ValueError(f"Unknown align mode: {align}")


def compute_phase_offsets(lcs: list[Lightcurve], multiband: bool, align_mode):
    phase_offsets = {}

    if multiband:
        offset = dispatch_phase_offsets(lcs, lcs[0].fopt.value, align_mode)
        for lc in lcs:
            phase_offsets[lc.band] = offset
    else:
        for lc in lcs:
            phase_offsets[lc.band] = dispatch_phase_offsets([lc], lc.fopt.value, align_mode)

    if align_mode in ["mean", "median"]:
        # choose reference band
        ref_band = lcs[0].band

        # compute reference slope
        ref_lc = next(lc for lc in lcs if lc.band == ref_band)

        phase_ref = ref_lc.phase
        phase_ref = (phase_ref - phase_offsets[ref_band]) % 1

        slope_ref = estimate_slope(phase_ref, ref_lc.brightness)

        for lc in lcs:
            band = lc.band
            phase = (lc.phase - phase_offsets[band]) % 1

            slope = estimate_slope(phase, lc.brightness)

            # flip if opposite orientation
            if slope * slope_ref < 0:
                phase_offsets[band] = (phase_offsets[band] + 0.5) % 1

    return phase_offsets


def handle_fold_arguments(lc: Lightcurve, **kwargs):
    lc = copy.deepcopy(lc)

    # align bands in phase
    aligned_phase = (lc.phase - kwargs["phase_offsets"][lc.band]) % 1
    setattr(lc, "phase", aligned_phase)

    # zero-align magnitudes
    subtract = kwargs.get("subtract", "median")
    if subtract == "mean":
        weights = 1 / np.power(lc.brightness_err, 2)
        brightness = lc.brightness - np.average(lc.brightness, weights=weights)
        setattr(lc, lc.brightness_type, brightness)
    elif subtract == "median":
        brightness = lc.brightness - np.median(lc.brightness)
        setattr(lc, lc.brightness_type, brightness)
    elif subtract is None:
        pass
    else:
        raise ValueError(f"Unexpected subtract mode '{subtract}'.")

    # repeat light curve
    phase, brightness, brightness_err = repeat_arrs(lc.phase, lc.brightness, lc.brightness_err, repeat=kwargs.get("repeat", 2))
    lc.phase = phase
    setattr(lc, lc.brightness_type, brightness)
    setattr(lc, f"{lc.brightness_type}_err", brightness_err)

    return lc
