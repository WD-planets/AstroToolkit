from bokeh.plotting import figure

from ...structures.definitions import Spectrum
from ..formatting import format_plot


def plot(spectrum: Spectrum, *args: any, **kwargs: any):
    """
    Plots an ATK Spectrum object
    """

    plot = figure(
        width=400,
        height=400,
        title=f"{spectrum.survey} Spectrum",
        x_axis_label="Wavelength / \u212b",
        y_axis_label=r"\[\text{flux / }10^{-17}\text{ erg}\text{cm}^{-2}\text{s}^{-1}\]" + "\u212b" + r"\[\:\:^{-1}\]",
        tools=("pan,wheel_zoom,box_zoom,reset"),
    )

    # plot spectrum
    plot.line(spectrum.wavelength, spectrum.flux, color="black", line_width=1)

    return format_plot("spectrum", plot)
