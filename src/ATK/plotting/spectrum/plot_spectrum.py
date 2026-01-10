from bokeh.plotting import figure

from ...structures.definitions import Spectrum
from ..formatting import format_plot
from ..plotting_core import get_axis_label


def plot(spectrum: Spectrum, *args: any, **kwargs: any):
    """
    Plots an ATK Spectrum object
    """

    plot = figure(
        width=400,
        height=400,
        title=f"{spectrum.survey} Spectrum",
        x_axis_label=get_axis_label(spectrum, "wavelength"),
        y_axis_label=get_axis_label(spectrum, "flux"),
        tools=("pan,wheel_zoom,box_zoom,reset"),
    )

    # plot spectrum
    plot.line(spectrum._get_attr_value("wavelength"), spectrum._get_attr_value("flux"), color="black", line_width=1)

    return format_plot("spectrum", plot)
