import warnings

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.stats import sigma_clip
from astropy.units import Quantity
from bokeh.io import show
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.signal import find_peaks

from ....plotting.formatting import format_plot
from ....structures.Spectrum import Spectrum

warnings.simplefilter("ignore", category=OptimizeWarning)

DEBUG = False

C_KMS = 299792.458
GAUSS_PARAMS = ("h", "a", "mu", "sigma")

INIT_WINDOW_WIDTH = 75  # initial window width in pixels
MIN_PROM = 1.0  # minimum prominence of peaks
MIN_SEP = 10  # minimum separation of 'new' peaks in pixels
MIN_WIDTH = 10

PROM_WEIGHT = 2.0  # weighting on prominence in scoring peaks
SEP_WEIGHT = 1.0  # weighting on separation in scoring peaks

INIT_LOCAL_WIN_WIDTH = 200  # width of initial local window in pixels
FINAL_LOCAL_WIN_WIDTH = 5  # width of final local window in sigma

CORE_SIGMA_WIDTH = 3  # width of 'core' of peak in sigma, used to determine edge pixels
MIN_EDGE_POINTS = 4  # minimum number of points to use in continuum estimation

FINAL_PLOTTING_WIDTH = 4  # width of final plotted Gaussian in sigma


def get_velocities(wav: Quantity, wav_ref: Quantity) -> np.ndarray:
    return (wav - wav_ref) / wav_ref * C_KMS * u.Unit("km s-1")


def gaussian(x, h, a, mu, sigma, sign=-1.0):
    return h + sign * a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def make_gaussian(fixed):
    def model(x, *free):
        params = {}
        i = 0
        for name in GAUSS_PARAMS:
            if name in fixed:
                params[name] = fixed[name]
            else:
                params[name] = free[i]
                i += 1
        return gaussian(x, **params)

    return model


def estimate_continuum(plot, x, y, mu, sigma):
    wing_mask = np.abs(x - mu) > CORE_SIGMA_WIDTH * sigma

    x_w = x[wing_mask]
    y_w = y[wing_mask]

    clipped = sigma_clip(y_w, sigma=3.0)

    if wing_mask.sum() < MIN_EDGE_POINTS:
        return np.median(y), x_w[~clipped.mask], y_w[~clipped.mask]

    if clipped.count() < MIN_EDGE_POINTS:
        return np.median(y_w), x_w[~clipped.mask], y_w[~clipped.mask]

    return np.median(clipped.data[~clipped.mask]), x_w[~clipped.mask], y_w[~clipped.mask]


def fit_gaussian(plot: figure, spectrum: Spectrum, p0: dict, snr_limit: float, fixed: dict | None = None):
    fixed = fixed or {}

    d_pix = np.mean(np.diff(spectrum.wavelength))

    init_pad = INIT_LOCAL_WIN_WIDTH * d_pix
    local_mask = (spectrum.wavelength >= p0["mu"] - init_pad) & (spectrum.wavelength <= p0["mu"] + init_pad)
    local_x, local_y = spectrum.wavelength[local_mask], spectrum.flux[local_mask]

    # if DEBUG:
    #     box = BoxAnnotation(left=p0["mu"] - INIT_PAD, right=p0["mu"] + INIT_PAD, fill_alpha=0.15, fill_color="orange")
    #     plot.add_layout(box)

    p0["sigma"] = 3 * d_pix
    p0["h"] = np.median(local_y)

    if p0["sign"] == -1:
        p0["a"] = np.abs(p0["h"] - local_y.min())
    else:
        p0["a"] = np.abs(local_y.max() - p0["h"])

    model = make_gaussian(fixed)
    free_p0 = [p0[k] for k in GAUSS_PARAMS if k not in fixed]

    popt, pcov = curve_fit(model, local_x, local_y, p0=free_p0, maxfev=20000)

    h, a, mu, sigma = popt

    sigma = abs(sigma)
    if not np.isfinite(sigma) or sigma < d_pix:
        sigma = FINAL_LOCAL_WIN_WIDTH * d_pix

    final_pad = FINAL_LOCAL_WIN_WIDTH * sigma
    local_mask = (spectrum.wavelength >= mu - final_pad) & (spectrum.wavelength <= mu + final_pad)
    local_x, local_y = spectrum.wavelength[local_mask], spectrum.flux[local_mask]

    # if DEBUG:
    #     box = BoxAnnotation(left=mu - FINAL_PAD, right=mu + FINAL_PAD, fill_alpha=0.15, fill_color="red")
    #     plot.add_layout(box)

    model_free = make_gaussian({})
    y_model = model_free(local_x, h, a, mu, sigma)
    y_resid = local_y - (y_model - h)

    continuum_h, wing_x, wing_y = estimate_continuum(plot, local_x, y_resid, mu, sigma)
    model = make_gaussian(fixed={"h": continuum_h})

    a0 = a
    mu0 = mu
    sigma0 = max(sigma, d_pix)

    popt, pcov = curve_fit(model, local_x, local_y, p0=[a0, mu0, sigma0], maxfev=20000)

    params = {"h": continuum_h}

    for name, val in zip(("a", "mu", "sigma"), popt):
        params[name] = val

    x_mask = (spectrum.wavelength >= mu - FINAL_PLOTTING_WIDTH * sigma) & (spectrum.wavelength <= mu + FINAL_PLOTTING_WIDTH * sigma)
    final_x = spectrum.wavelength[x_mask]
    final_y = model(final_x, *popt)

    local_y_smooth = gaussian_filter1d(local_y, sigma=5)
    resid = local_y - local_y_smooth
    noise = 1.4826 * np.median(np.abs(resid - np.median(resid)))

    peak_amp = final_y.max() - final_y.min()
    peak_width = final_x.max() - final_x.min()
    peak_pix = len(final_x)

    peak_snr = peak_amp / noise

    filtering = {
        "peak_width": [peak_width],
        "peak_amp": peak_amp,
        "peak_snr": peak_snr,
        "peak_pix": peak_pix,
        "local_noise": noise,
        "local_std": np.std(final_y),
        "continuum_h": continuum_h,
    }

    if peak_width > 0.1 * (spectrum.wavelength.max() - spectrum.wavelength.min()):
        return None, f"width too wide ({peak_width:.2f})", filtering

    if DEBUG:
        plot.line(local_x, [continuum_h] * len(local_x), line_color="red", legend_label="Continuum Level")
        plot.scatter(wing_x, wing_y, size=4, alpha=0.6, color="orange", legend_label="Continuum Points")

    if peak_width < 10:
        return None, f"width too narrow ({peak_pix} px)", filtering

    if peak_snr < snr_limit:
        return None, f"snr too low ({peak_snr:.2f})", filtering

    final_model = {
        "mu": params["mu"],
        "sigma": params["sigma"],
        "a": params["a"],
        "h": params["h"],
        "sign": p0["sign"],
        "model_x": final_x,
        "model_y": final_y,
    }

    return final_model, "accepted", filtering


def detect_features(plot, spectrum: Spectrum, min_prominence: float, smoothing: int):
    flux_smooth = gaussian_filter1d(spectrum.flux, sigma=smoothing)

    # emission peaks
    em_peaks, em_props = find_peaks(flux_smooth, prominence=min_prominence, width=MIN_WIDTH)

    # absorption peaks
    abs_peaks, abs_props = find_peaks(-flux_smooth, prominence=min_prominence, width=MIN_WIDTH)

    # combine
    peaks = np.concatenate([em_peaks, abs_peaks])
    signs = np.concatenate([np.ones(len(em_peaks), dtype=int), -np.ones(len(abs_peaks), dtype=int)])
    prominences = np.concatenate([em_props["prominences"], abs_props["prominences"]])
    widths = np.concatenate([em_props["widths"], abs_props["widths"]])

    score = PROM_WEIGHT * prominences + SEP_WEIGHT * widths

    data = pd.DataFrame({"idx": peaks, "sign": signs, "prominence": prominences, "width": widths, "score": score}).sort_values("idx")

    keep = []
    last_idx = None
    for _, row in data.iterrows():
        if last_idx is None:
            keep.append(row)
            last_idx = row["idx"]
            continue

        if row["idx"] - last_idx < MIN_SEP:
            if row["score"] > keep[-1]["score"]:
                keep[-1] = row
        else:
            keep.append(row)
            last_idx = row["idx"]

    data = pd.DataFrame(keep)

    data["vel"] = spectrum.wavelength[data["idx"].to_numpy(dtype=int)]
    data["flux"] = spectrum.flux[data["idx"].to_numpy(dtype=int)]

    return flux_smooth, data


def subtract_peak(wavelength: np.ndarray, flux: np.ndarray, model_x: np.ndarray, model_y: np.ndarray, h: float):
    mask = (wavelength >= model_x.min()) & (wavelength <= model_x.max())
    flux_new = flux.copy()
    flux_new[mask] -= model_y - h

    return flux_new


def mad_sigma(x):
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def remove_false_features(plot, spectrum, fits):
    from scipy.ndimage import generic_filter, median_filter

    resid_flux = spectrum.flux.copy()
    for fit in fits:
        resid_flux = subtract_peak(spectrum.wavelength, resid_flux, fit["model_x"], fit["model_y"], fit["h"])

    resid_flux = sigma_clip(resid_flux, sigma=3)
    resid_flux = resid_flux.filled(np.nan)
    resid_flux = np.asarray(resid_flux, dtype=float)

    window = 200
    continuum = median_filter(resid_flux, size=window)

    sigma = generic_filter(resid_flux, mad_sigma, size=window)

    keep = []
    for fit in fits:
        idx = np.abs(spectrum.wavelength - fit["mu"]).argmin()
        peak_val = fit["model_y"].max() if fit["sign"] == 1 else fit["model_y"].min()

        if np.abs(peak_val - continuum[idx]) > 3 * sigma[idx]:
            keep.append(fit)

    if DEBUG:
        plot.line(spectrum.wavelength, resid_flux, line_color="blue", legend_label="Residual Flux")
        plot.line(spectrum.wavelength, continuum, line_color="red", legend_label="Running Median Flux")

    return keep


def do_fitting(spectrum: Spectrum, prominence: float, snr: float, smoothing: int, debug=False):
    global DEBUG
    DEBUG = debug

    wavelength_unit = spectrum.wavelength.unit

    spectrum.wavelength = spectrum.wavelength.value
    spectrum.flux = spectrum.flux.value

    plot = figure(width=1000, height=500, x_axis_label=f"Wavelength / {wavelength_unit.to_string('unicode')}", y_axis_label="Flux")
    plot.line(spectrum.wavelength, spectrum.flux, line_color="black", line_alpha=0.5)

    flux_smooth, peak_data = detect_features(plot, spectrum, prominence, smoothing)

    fits, responses = [], []
    filtering = []
    for index, row in peak_data.iterrows():
        p0 = {}
        p0["mu"] = row["vel"]
        p0["sign"] = row["sign"]

        fit, response, filter = fit_gaussian(plot, spectrum, p0, snr_limit=snr, fixed={})

        if fit:
            fits.append(fit)
        responses.append(response)
        filtering.append(filter)

    fits = remove_false_features(plot, spectrum, fits)
    for fit in fits:
        plot.line(fit["model_x"], fit["model_y"], line_color="limegreen", line_width=2, legend_label="Detected Peaks")

    plot.line(spectrum.wavelength, flux_smooth, line_alpha=0.3, line_width=3, line_color="red", legend_label="Smoothed Flux")

    filtering_df = pd.concat(pd.DataFrame(dct) for dct in filtering).reset_index(drop=True)

    if DEBUG:
        df = pd.DataFrame(
            {
                "x": spectrum.wavelength[peak_data["idx"].to_numpy(dtype=int)],
                "y": flux_smooth[peak_data["idx"].to_numpy(dtype=int)],
                "response": responses,
            }
        )
        df = pd.concat([df, filtering_df], axis=1)
        source = ColumnDataSource(df)
        scatter = plot.scatter("x", "y", source=source, marker="x", size=10, color="red", legend_label="Detected Peaks")
        hvr = HoverTool(tooltips=[(name, f"@{name}") for name in df.columns.values.tolist()])
        hvr.renderers = [scatter]
        plot.add_tools(hvr)

    plot = format_plot("spectrum", plot)

    show(plot)
