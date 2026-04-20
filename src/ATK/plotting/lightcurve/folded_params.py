import copy

import astropy.units as u
import numpy as np


def repeat_arrs(x, *ys, repeat=2):
    x_unit = getattr(x, "unit", None)
    y_units = [getattr(y, "unit", None) for y in ys]

    x_val = x.to_value() if x_unit else x
    ys_val = [y.to_value() if getattr(y, "unit", None) else y for y in ys]

    xs = []
    ys_out = [[] for _ in ys]

    for k in range(repeat):
        xs.append(x_val + k)
        for i, y in enumerate(ys_val):
            ys_out[i].append(y)

    x_rep = np.concatenate(xs)
    ys_rep = [np.concatenate(ylist) for ylist in ys_out]

    # reattach units
    if x_unit:
        x_rep = x_rep * x_unit
    ys_rep = [y * unit if unit else y for y, unit in zip(ys_rep, y_units)]

    return (x_rep, *ys_rep)


def compute_band_offsets(lcs, align_method="median"):
    ref_lc = lcs[0]

    if align_method is None:
        return {lc.band: 0.0 * ref_lc._brightness.unit for lc in lcs}

    if align_method not in ["mean", "median"]:
        raise ValueError(f"Unexpected align method '{align_method}'.")

    if align_method == "mean":
        ref_val = np.average(ref_lc._brightness, weights=1 / ref_lc._brightness_err**2)
    else:
        ref_val = np.median(ref_lc._brightness)

    band_offsets = {}

    for lc in lcs:
        if align_method == "mean":
            val = np.average(lc._brightness, weights=1 / lc._brightness_err**2)
        else:
            val = np.median(lc._brightness)

        band_offsets[lc.band] = ref_val - val  # keeps units

    return band_offsets


def compute_phase_offset_midpoint(lcs, freq, align_method="median", nbins=50):
    phase = np.concatenate([lc.phase.to_value(1) for lc in lcs])
    flux = np.concatenate([lc._brightness.value for lc in lcs])

    midpoint = np.mean(flux) if align_method == "mean" else np.median(flux)

    bins = np.linspace(0, 1, nbins + 1)
    digitized = np.digitize(phase, bins)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.full(nbins, np.nan)

    for i in range(1, nbins + 1):
        vals = flux[digitized == i]
        if len(vals):
            bin_means[i - 1] = np.nanmean(vals)

    valid = ~np.isnan(bin_means)
    if not np.any(valid):
        return 0.0

    idx = np.where(valid)[0][np.argmin(np.abs(bin_means[valid] - midpoint))]

    return float(bin_centers[idx])


def compute_phase_offset_extrema(lcs, freq, mode="max"):
    brightness_type = list(set(lc._brightness_type for lc in lcs))
    if len(brightness_type) > 1:
        raise ValueError("Mixed brightness types.")
    brightness_type = brightness_type[0]

    phase = np.concatenate([lc.phase.to_value(1) for lc in lcs])
    flux = np.concatenate([lc._brightness.value for lc in lcs])

    bins = np.linspace(0, 1, 50)
    digitized = np.digitize(phase, bins)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_means = np.array([np.nanmean(flux[digitized == i]) if np.any(digitized == i) else np.nan for i in range(1, len(bins))])

    valid = ~np.isnan(bin_means)
    phase = bin_centers[valid]
    flux = bin_means[valid]

    idx = np.argsort(phase)
    phase, flux = phase[idx], flux[idx]

    phase_ext = np.concatenate([phase - 1, phase, phase + 1])
    flux_ext = np.tile(flux, 3)

    idx = np.argsort(phase_ext)
    phase_ext, flux_ext = phase_ext[idx], flux_ext[idx]

    window = max(5, len(flux) // 20)
    kernel = np.ones(window) / window
    smooth = np.convolve(flux_ext, kernel, mode="same")

    mask = (phase_ext >= 0) & (phase_ext < 1)
    phase, smooth = phase_ext[mask], smooth[mask]

    if mode == "max":
        i = np.argmin(smooth) if brightness_type == "mag" else np.argmax(smooth)
    elif mode == "min":
        i = np.argmax(smooth) if brightness_type == "mag" else np.argmin(smooth)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return float(phase[i])


def estimate_slope(phase, flux, window=0.1):
    flux = flux.value

    phase_ext = np.concatenate([phase - 1, phase, phase + 1])
    flux_ext = np.tile(flux, 3)

    mask = (phase_ext >= -window) & (phase_ext <= window)

    if np.sum(mask) < 5:
        return 0.0

    coeffs = np.polyfit(phase_ext[mask], flux_ext[mask], 1)

    return float(coeffs[0])


def dispatch_phase_offsets(lcs, freq, align="median"):
    if align in ["mean", "median"]:
        return compute_phase_offset_midpoint(lcs, freq, align)
    elif align in ["max", "min"]:
        return compute_phase_offset_extrema(lcs, freq, align)
    else:
        raise ValueError(f"Unknown align mode: {align}")


def compute_phase_offsets(lcs, multiband, align_mode):
    phase_offsets = {}

    if multiband:
        offset = dispatch_phase_offsets(lcs, lcs[0].fopt, align_mode)
        for lc in lcs:
            phase_offsets[lc.band] = offset
    else:
        for lc in lcs:
            phase_offsets[lc.band] = dispatch_phase_offsets([lc], lc.fopt, align_mode)

    if align_mode in ["mean", "median"]:
        ref_band = lcs[0].band
        ref_lc = next(lc for lc in lcs if lc.band == ref_band)

        phase_ref = (ref_lc.phase.to_value(1) - phase_offsets[ref_band]) % 1
        slope_ref = estimate_slope(phase_ref, ref_lc._brightness)

        for lc in lcs:
            phase = (lc.phase.to_value(1) - phase_offsets[lc.band]) % 1
            slope = estimate_slope(phase, lc._brightness)

            if slope * slope_ref < 0:
                phase_offsets[lc.band] = (phase_offsets[lc.band] + 0.5) % 1

    return phase_offsets


def handle_fold_arguments(lc, **kwargs):
    lc = copy.deepcopy(lc)

    offset = kwargs["phase_offsets"][lc.band]

    phase = lc.phase - offset * u.one
    phase = (phase.to_value(u.one) % 1) * u.one

    flux = lc._brightness
    err = lc._brightness_err

    subtract = kwargs.get("subtract", "median")

    if subtract == "mean":
        weights = 1 / err**2
        flux = flux - np.average(flux, weights=weights)
    elif subtract == "median":
        flux = flux - np.median(flux)
    elif subtract is not None:
        raise ValueError(f"Unexpected subtract mode '{subtract}'.")

    phase, flux, err = repeat_arrs(phase, flux, err, repeat=kwargs.get("repeat", 2))

    lc.phase = phase
    setattr(lc, lc._brightness_type, flux)
    setattr(lc, f"{lc._brightness_type}_err", err)

    return lc
