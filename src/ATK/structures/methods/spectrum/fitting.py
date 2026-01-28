import numpy as np
from astropy.units import Quantity
from bokeh.io import show
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from scipy.signal import find_peaks

from ....structures.Spectrum import Spectrum

BINS = 1000


def get_velocities(wav: Quantity, wav_ref: Quantity) -> np.ndarray:
    return ((wav - wav_ref) / wav_ref) * 3e5


def basic_bin(data, width):
    return data[: (data.size // width) * width].reshape(-1, width).mean(axis=1)


def get_priors(spectrum: Spectrum, vel: np.ndarray, width: float):
    # bin spectrum to smooth out noise (need to figure out how to adapt this to amount of data/peak width to not smooth out a weak/thin peak)
    flux_b = basic_bin(spectrum.flux, BINS)
    vel_b = basic_bin(vel, BINS)

    # invert for absorption features
    vel_b = -vel_b

    # find peaks
    peaks, _ = find_peaks(vel_b, distance=1)

    min_fluxes = flux_b[peaks]
    min_vel = vel_b[peaks]


def do_fitting(spectrum: Spectrum, wavelengths: list[float], widths: list[float]):
    # get velocities of requested features
    if not isinstance(wavelengths, np.ndarray):
        wavelengths = np.asarray(wavelengths)

    line_vels = get_velocities(wavelengths, wavelengths[0])
    vel = get_velocities(spectrum.wavelength, wavelengths[0])

    fig = figure(
        width=1000,
        height=500,
        x_axis_label=r"Velocity / $$\text{km\,s}^{-1}$$",
        y_axis_label=r"Flux / $$10^{-16}\text{erg}\,\text{s}^{-1}\,\text{cm}^{-2}\,\text{Angstrom}^{-1}$$",
    )
    data = ColumnDataSource(data={"v": vel, "f": spectrum.flux})
    fig.line(x="v", y="f", source=data, line_alpha=0.5)

    show(fig)
