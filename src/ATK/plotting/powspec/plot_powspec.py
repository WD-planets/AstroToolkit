import numpy as np
from bokeh.models import Range1d
from bokeh.plotting import figure

from ...structures.Powspec import Powspec
from ..formatting import format_plot


def plot(pspec: Powspec, *args: tuple, **kwargs: dict):
    plot = figure(
        width=400,
        height=400,
        x_axis_label=r"Frequency / \[\text{days}^{-1}\]",
        y_axis_label="Lomb-Scargle Power",
        title=f"{pspec.survey} {pspec.band} L-S Power Spectrum",
        tools=("pan,wheel_zoom,box_zoom,reset"),
    )

    plot.line(x=pspec.frequency, y=pspec.power)

    plot.y_range = Range1d(0, np.nanmax(pspec.power) * 1.1)
    plot.x_range = Range1d(np.nanmin(pspec.frequency.value), np.nanmax(pspec.frequency.value))

    pspec._plot_id = plot.id

    return format_plot("powspec", plot)
