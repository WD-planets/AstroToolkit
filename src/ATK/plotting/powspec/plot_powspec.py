import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, Range1d
from bokeh.plotting import figure

from ...structures.Powspec import Powspec
from ..formatting import format_plot
from ..plotting_core import get_axis_label


def plot(pspec: Powspec, *args: tuple, **kwargs: dict):
    plot = figure(
        width=400,
        height=400,
        x_axis_label=get_axis_label(pspec, "frequency"),
        y_axis_label="Lomb-Scargle Power",
        title=f"{pspec.survey} {pspec.band} L-S Power Spectrum",
        tools=("pan,wheel_zoom,box_zoom,reset"),
    )

    hvr = HoverTool(tooltips=[("fopt", f"@fopt {pspec.fopt.unit}"), ("popt", f"@popt {pspec.popt.unit}")])

    df = pd.DataFrame({"freq": pspec.frequency, "power": pspec.power, "fopt": pspec.fopt, "popt": pspec.popt})
    source = ColumnDataSource(df)
    line = plot.line(x="freq", y="power", source=source, legend_label=f"{pspec.band}-band Power")
    hvr.renderers = [line]
    plot.add_tools(hvr)

    plot.y_range = Range1d(0, np.nanmax(pspec.power) * 1.1)
    plot.x_range = Range1d(np.nanmin(pspec.frequency.value), np.nanmax(pspec.frequency.value))

    pspec._plot_id = plot.id

    return format_plot("powspec", plot)
