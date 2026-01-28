import numpy as np
from bokeh.io import show
from bokeh.plotting import figure
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from ....structures.Spectrum import Spectrum

C_KMS = 299792.458  # speed of light in km/s
GAUSS_PARAMS = ("H", "A", "mu", "sigma")
BINS = 1000


def velocity_from_wavelength(wav, wav_ref):
    return (wav - wav_ref) / wav_ref * C_KMS


def gaussian(x, H, A, mu, sigma, sign=-1.0):
    """
    General Gaussian. sign=-1 absorption, +1 emission.
    """

    return H + sign * A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def make_gaussian(fixed, sign):
    def model(x, *free):
        params = {}
        i = 0
        for name in GAUSS_PARAMS:
            if name in fixed:
                params[name] = fixed[name]
            else:
                params[name] = free[i]
                i += 1
        return gaussian(x, **params, sign=sign)

    return model


def fit_gaussian(x, y, p0, fixed=None, absorption=True):
    fixed = fixed or {}
    sign = -1.0 if absorption else 1.0

    model = make_gaussian(fixed, sign)
    free_p0 = [p0[k] for k in GAUSS_PARAMS if k not in fixed]

    popt, pcov = curve_fit(model, x, y, p0=free_p0, maxfev=20_000)

    params = fixed.copy()
    i = 0
    for name in GAUSS_PARAMS:
        if name not in params:
            params[name] = popt[i]
            i += 1

    return params, model(x, *popt)


def subset(x, y, lo, hi):
    mask = (x >= lo) & (x <= hi)
    return x[mask], y[mask]


def bin_mean(arr, width):
    n = arr.size // width
    return arr[: n * width].reshape(n, width).mean(axis=1)


def estimate_priors(vel, flux, feature_width, bins=BINS):
    vel_b = bin_mean(vel, bins)
    flux_b = bin_mean(flux, bins)

    # absorption → minima
    peaks, _ = find_peaks(flux_b, prominence=np.std(flux_b))

    if len(peaks) == 0:
        raise RuntimeError("No absorption feature found")

    idx = peaks[np.argmin(flux_b[peaks])]
    mu = vel_b[idx]
    A = np.mean(flux_b) - flux_b[idx]
    H = np.median(flux_b)

    # convert FWHM-like width → sigma
    sigma = feature_width / (2 * np.sqrt(2 * np.log(2)))

    return dict(H=H, A=A, mu=mu, sigma=sigma)


def do_fitting(spectrum: Spectrum, wavelengths, widths):
    wavelengths = np.asarray(wavelengths)
    vel = velocity_from_wavelength(spectrum.wavelength, wavelengths[0])

    fig = figure(width=1000, height=500, x_axis_label=r"Velocity / km s$^{-1}$", y_axis_label=r"Flux")
    fig.line(vel, spectrum.flux, alpha=0.5)

    results = []

    for wav, width in zip(wavelengths, widths):
        priors = estimate_priors(vel, spectrum.flux, width)

        lo, hi = priors["mu"] - width, priors["mu"] + width
        xsub, ysub = subset(vel, spectrum.flux, lo, hi)

        params, fit_y = fit_gaussian(xsub, ysub, priors, absorption=True)

        print(xsub, fit_y)

        fig.line(xsub, fit_y, line_color="green", line_width=2)
        results.append(params)

    show(fig)
    return results
