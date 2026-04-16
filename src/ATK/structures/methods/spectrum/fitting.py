from __future__ import annotations

import copy
import warnings
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from astropy.stats import sigma_clip
from astropy.units import Quantity
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from scipy.ndimage import gaussian_filter1d, generic_filter, median_filter
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.signal import find_peaks

from ....plotting.formatting import format_plot

if TYPE_CHECKING:
    from ....structures.Spectrum import Spectrum

warnings.simplefilter("ignore", category=OptimizeWarning)

C_KMS = 299792.458
GAUSS_PARAMS = ("h", "a", "mu", "sigma")

MIN_SEP = 10  # minimum separation below which peaks are discarded
MIN_WIDTH = 10  # minimum width of a peak in pixels

PROM_WEIGHT = 2.0  # weighting on prominence in scoring peaks
SEP_WEIGHT = 1.0  # weighting on separation in scoring peaks

INIT_LOCAL_WIN_WIDTH = 200  # width of initial local window in wavelength
FINAL_LOCAL_WIN_WIDTH = 5  # width of final local window in sigma

CORE_SIGMA_WIDTH = 3  # width of 'core' of peak in sigma, used to determine edge pixels
MIN_EDGE_POINTS = 4  # minimum number of points to use in continuum estimation

FINAL_PLOTTING_WIDTH = 4  # width of final plotted Gaussian in sigma


def gaussian(x: np.ndarray, h: float, a: float, mu: float, sigma: float, sign: int = -1):
    return h + sign * a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def make_gaussian(fixed):
    """
    Models a Gaussian with or without fixed parameters (h/a/mu/sigma)
    """

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


def estimate_continuum(plot: figure, x: np.ndarray, y: np.ndarray, mu: float, sigma: float):
    """
    Estimates the local continuum by masking out a feature
    """

    wing_mask = np.abs(x - mu) > CORE_SIGMA_WIDTH * sigma

    x_w = x[wing_mask]
    y_w = y[wing_mask]

    if wing_mask.sum() < MIN_EDGE_POINTS:
        return np.median(y), x_w, y_w

    clipped = sigma_clip(y_w, sigma=3.0)

    if clipped.count() < MIN_EDGE_POINTS:
        return np.median(y_w), x_w[~clipped.mask], y_w[~clipped.mask]

    return np.median(clipped.data[~clipped.mask]), x_w[~clipped.mask], y_w[~clipped.mask]


def fit_gaussian(plot: figure, spectrum: Spectrum, p0: dict, snr_limit: float):
    """
    Fits Gaussians to spectral features in two stages
    """

    # get width of one sample (pixel)
    d_pix = np.mean(np.diff(spectrum.wavelength))

    # define initial locoal window
    init_pad = INIT_LOCAL_WIN_WIDTH * d_pix
    local_mask = (spectrum.wavelength >= p0["mu"] - init_pad) & (spectrum.wavelength <= p0["mu"] + init_pad)
    local_x, local_y = spectrum.wavelength[local_mask], spectrum.flux[local_mask]

    # set up priors
    p0["sigma"] = 3 * d_pix
    p0["h"] = np.median(local_y)
    p0["a"] = np.abs(p0["h"] - local_y.min()) if p0["sign"] == -1 else np.abs(local_y.max() - p0["h"])

    # set up basic Gaussian model
    model = make_gaussian({})
    free_p0 = [p0[k] for k in GAUSS_PARAMS]

    # get initial fit parameters
    popt, pcov = curve_fit(model, local_x, local_y, p0=free_p0, maxfev=20000)
    h, a, mu, sigma = popt

    # get peak width
    sigma = abs(sigma)
    if not np.isfinite(sigma) or sigma < d_pix:
        sigma = FINAL_LOCAL_WIN_WIDTH * d_pix

    # set up secondary window as multiple of peak width
    final_pad = FINAL_LOCAL_WIN_WIDTH * sigma
    local_mask = (spectrum.wavelength >= mu - final_pad) & (spectrum.wavelength <= mu + final_pad)
    local_x, local_y = spectrum.wavelength[local_mask], spectrum.flux[local_mask]

    # create model from fit parameters and subtract from local flux
    model_free = make_gaussian({})
    y_model = model_free(local_x, h, a, mu, sigma)
    y_resid = local_y - (y_model - h)

    # estimate continuum and get a new model with fixed height
    continuum_h, wing_x, wing_y = estimate_continuum(plot, local_x, y_resid, mu, sigma)
    model = make_gaussian(fixed={"h": continuum_h})

    # set up other priors
    a0 = a
    mu0 = mu
    sigma0 = max(sigma, d_pix)

    # get final fit
    popt, pcov = curve_fit(model, local_x, local_y, p0=[a0, mu0, sigma0], maxfev=20000)

    # set up params dict
    params = {"h": continuum_h}
    for name, val in zip(("a", "mu", "sigma"), popt):
        params[name] = val

    # clip final x for plotting and get final y
    x_mask = (spectrum.wavelength >= mu - FINAL_PLOTTING_WIDTH * sigma) & (spectrum.wavelength <= mu + FINAL_PLOTTING_WIDTH * sigma)
    final_x = spectrum.wavelength[x_mask]
    final_y = model(final_x, *popt)

    # smooth local y and subtract, use MAD to estimate local noise
    local_y_smooth = gaussian_filter1d(local_y, sigma=5)
    resid = local_y - local_y_smooth
    noise = 1.4826 * np.median(np.abs(resid - np.median(resid)))

    # get filtering parameters
    peak_amp = final_y.max() - final_y.min()
    peak_width = final_x.max() - final_x.min()
    peak_pix = len(final_x)
    peak_snr = peak_amp / noise

    # used for HoverTool when DEBUG is True
    filtering = {
        "peak_width": [peak_width],
        "peak_amp": peak_amp,
        "peak_snr": peak_snr,
        "peak_pix": peak_pix,
        "local_noise": noise,
        "local_std": np.std(final_y),
        "continuum_h": continuum_h,
    }

    # discard peaks that cover far too much of the spectrum
    if peak_width > 0.1 * (spectrum.wavelength.max() - spectrum.wavelength.min()):
        return None, f"width too wide ({peak_width:.2f})", filtering

    if DEBUG:
        plot.line(local_x, [continuum_h] * len(local_x), line_color="red", legend_label="Continuum Level")
        plot.scatter(wing_x, wing_y, size=4, alpha=0.6, color="orange", legend_label="Continuum Points")

    # discard really thin peaks
    if peak_width < 10:
        return None, f"width too narrow ({peak_pix} px)", filtering

    # discard peaks with poor snr
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
        "snr": peak_snr,
    }

    return final_model, "accepted", filtering


def detect_features(plot: figure, spectrum: Spectrum, min_prominence: float, smoothing: int):
    """
    Identifies features in a spectrum
    """

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

    # generate score from prominence and width
    score = PROM_WEIGHT * prominences + SEP_WEIGHT * widths

    data = pd.DataFrame({"idx": peaks, "sign": signs, "prominence": prominences, "width": widths, "score": score}).sort_values("idx")

    keep = []
    last_idx = None
    for _, row in data.iterrows():
        if last_idx is None:
            keep.append(row)
            last_idx = row["idx"]
            continue

        # keep only best of two close peaks
        if row["idx"] - last_idx < MIN_SEP:
            if row["score"] > keep[-1]["score"]:
                keep[-1] = row
        else:
            keep.append(row)
            last_idx = row["idx"]

    data = pd.DataFrame(keep)

    data["vel"] = spectrum.wavelength[data["idx"].to_numpy(dtype=int)]
    data["flux"] = spectrum.flux[data["idx"].to_numpy(dtype=int)]

    return data


def subtract_peak(wavelength: np.ndarray, flux: np.ndarray, model_x: np.ndarray, model_y: np.ndarray, h: float):
    """
    Subtracts a peak from a spectrum given its height
    """

    mask = (wavelength >= model_x.min()) & (wavelength <= model_x.max())
    flux_new = flux.copy()
    flux_new[mask] -= model_y - h

    return flux_new


def mad_sigma(x):
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def remove_false_features(plot, spectrum, fits):
    """
    Removes false features created by the intersection of two close real features
    """

    # subtract features from spectrum
    resid_flux = spectrum.flux.copy()
    for fit in fits:
        resid_flux = subtract_peak(spectrum.wavelength, resid_flux, fit["model_x"], fit["model_y"], fit["h"])

    # sigma clip
    resid_flux = sigma_clip(resid_flux, sigma=3)
    resid_flux = resid_flux.filled(np.nan)
    resid_flux = np.asarray(resid_flux, dtype=float)

    # running median window
    window = 200
    continuum = median_filter(resid_flux, size=window)

    # running Median Absolute Deviation window to get width of continuum
    sigma = generic_filter(resid_flux, mad_sigma, size=window)

    keep = []
    for fit in fits:
        idx = np.abs(spectrum.wavelength - fit["mu"]).argmin()

        # get value at peak of feature
        peak_val = fit["model_y"].max() if fit["sign"] == 1 else fit["model_y"].min()

        # filter out features where the peak value is close to the continuum
        if np.abs(peak_val - continuum[idx]) > 3 * sigma[idx]:
            keep.append(fit)

    if DEBUG:
        plot.line(spectrum.wavelength, resid_flux, line_color="blue", legend_label="Residual Flux")
        plot.line(spectrum.wavelength, continuum, line_color="red", legend_label="Running Median Flux")

    return keep


def do_fitting(plot: figure, spectrum: Spectrum, prominence: float = 2, smoothing: int = 3, snr: float = 3, **kwargs) -> figure | Quantity:
    """
    Identifies and fits any number of spectral absorption and emission features
    """

    # set debugging mode
    global DEBUG
    DEBUG = kwargs.get("debug", False)

    # convert quantity arrays to basic arrays
    spectrum_copy = copy.deepcopy(spectrum)
    wavelength_unit = spectrum_copy.wavelength.unit
    spectrum_copy.wavelength = spectrum.wavelength.copy().value
    spectrum_copy.flux = spectrum.flux.copy().value

    if plot is None:
        plot = figure(width=1000, height=500, x_axis_label=f"Wavelength / {wavelength_unit.to_string('unicode')}", y_axis_label="Flux")

    # get features
    peak_data = detect_features(plot, spectrum_copy, prominence, smoothing)

    flux_smooth = gaussian_filter1d(spectrum_copy.flux, sigma=smoothing)

    # loop through returned peaks
    fits, responses, filtering = [], [], []
    for index, row in peak_data.iterrows():
        # assemble priors
        p0 = {}
        p0["mu"] = row["vel"]
        p0["sign"] = row["sign"]

        # fit gaussians to features
        fit, response, filter = fit_gaussian(plot, spectrum_copy, p0, snr_limit=snr)

        if fit:
            fits.append(fit)

        responses.append(response)
        filtering.append(filter)

    # remove fake features that result from two close real features
    fits = remove_false_features(plot, spectrum_copy, fits)

    # plot peaks
    features, peaks = [], []
    for fit in fits:
        plot.line(fit["model_x"], fit["model_y"], line_color="limegreen", line_width=2, legend_label="Detected Peaks")
        features.append(fit["mu"])
        peak_val = fit["model_y"].max() if fit["sign"] == 1 else fit["model_y"].min()
        peaks.append(peak_val)

    if kwargs.get("get_features"):
        return np.asarray(features), peaks

    # plot smoothed flux
    plot.line(spectrum_copy.wavelength, flux_smooth, line_alpha=0.3, line_width=3, line_color="red", legend_label="Smoothed Flux")

    # peak filtering overlay
    filtering_df = pd.concat(pd.DataFrame(dct) for dct in filtering).reset_index(drop=True)
    if DEBUG:
        df = pd.DataFrame(
            {
                "x": spectrum_copy.wavelength[peak_data["idx"].to_numpy(dtype=int)],
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

    return plot
