import numpy as np
import pandas as pd
from astropy.units import Quantity
from bokeh.io import show
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from ....structures.Spectrum import Spectrum

C_KMS = 299792.458
BIN_N = 10
GAUSS_PARAMS = ("h", "a", "mu", "sigma")


def get_velocities(wav: Quantity, wav_ref: Quantity) -> np.ndarray:
    return (wav - wav_ref) / wav_ref * C_KMS


def basic_bin(data: np.ndarray, bin_n: int):
    return data[: (data.size // bin_n) * bin_n].reshape(-1, bin_n).mean(axis=1)


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


def fit_gaussian(x, y, p0, fixed=None):
    fixed = fixed or {}

    model = make_gaussian(fixed)
    free_p0 = [p0[k] for k in GAUSS_PARAMS if k not in fixed]

    popt, pcov = curve_fit(model, x, y, p0=free_p0, maxfev=20000)

    params = fixed.copy()
    i = 0
    for name in GAUSS_PARAMS:
        if name not in params:
            params[name] = popt[i]
            i += 1

    return params, model(x, *popt)


def estimate_priors(plot: figure, spectrum: Spectrum, line_vel: float, width: float, completed: list[int]):
    flux_b = basic_bin(spectrum.flux, BIN_N)
    vel_b = basic_bin(spectrum.velocity, BIN_N)

    plot.line(vel_b, flux_b, line_alpha=0.25, line_width=5, line_color="black")

    # absorption → invert flux
    peaks, props = find_peaks(-flux_b, prominence=0.5 * np.std(flux_b))

    if len(peaks) == 0:
        return None

    hvr = HoverTool(tooltips=[("peak", "@peak"), ("prominence", "@prominence")])
    source = pd.DataFrame(
        {
            "peak": peaks,
            "prominence": props["prominences"],
            "vel": vel_b[peaks],
            "flux": flux_b[peaks],
            "sep": np.abs(vel_b[peaks] - line_vel),
        }
    )
    source = source.sort_values(by=["sep", "prominence"])
    scatter = plot.scatter(source=source, x="vel", y="flux", marker="x", size=10, line_color="red")
    hvr.renderers = [scatter]
    plot.add_tools(hvr)
    peaks = source["peak"].to_numpy()

    idx = peaks[0]
    while idx in completed:
        idx = peaks[0]
        peaks = [peak for peak in peaks if peak != idx]

    mu = vel_b[idx]
    h = np.median(flux_b)
    a = h - flux_b[idx]

    sigma = width / (2 * np.sqrt(2 * np.log(2)))

    return {"id": idx, "h": h, "a": a, "mu": mu, "sigma": sigma}


def do_fitting(spectrum: Spectrum, wavelengths: list[float], widths: list[float]):
    # get velocities of requested features
    if not isinstance(wavelengths, np.ndarray):
        wavelengths = np.asarray(wavelengths)

    line_vels = get_velocities(wavelengths, wavelengths[0])
    vel = get_velocities(spectrum.wavelength, wavelengths[0])

    vel_spec = spectrum.vspec(wav_ref=wavelengths[0], inplace=False)

    plot = figure(
        width=1000,
        height=500,
        x_axis_label=r"Velocity / $$\text{km\,s}^{-1}$$",
        y_axis_label=r"Flux / $$10^{-16}\text{erg}\,\text{s}^{-1}\,\text{cm}^{-2}\,\text{Angstrom}^{-1}$$",
    )

    data = ColumnDataSource(data={"v": vel, "f": spectrum.flux})
    plot.line(x="v", y="f", source=data, line_color="black")

    completed = []
    for line_vel, width in zip(line_vels, widths):
        priors = estimate_priors(plot, vel_spec, line_vel, width, completed)

        if priors is None:
            continue

        completed.append(priors["id"])

        low, high = priors["mu"] - width, priors["mu"] + width
        peak_data = vel_spec.crop(min=low, max=high, inplace=False)

        params, fit_y = fit_gaussian(peak_data.velocity, peak_data.flux, priors)
        plot.line(peak_data.velocity, fit_y, line_color="green", line_width=2)

    show(plot)
