import numpy as np
from bokeh.models import ColumnDataSource
from bokeh.models.formatters import BasicTickFormatter
from bokeh.plotting import figure
from bokeh.transform import linear_cmap

from ...structures.definitions import Lightcurve
from ..colours import assign_gradient_palettes
from ..formatting import format_plot


def plot_band(plot: figure, lc: Lightcurve, palette: list[str], time_min: float, time_format: str):
    # time handling
    time = lc.mjd
    if time_format == "reduced":
        time = [t - time_min for t in time]

    # get brightness and error columns
    y = getattr(lc, lc.brightness_type)
    y_err = getattr(lc, f"{lc.brightness_type}_err")

    source = ColumnDataSource(data={"time": time, "y": y})

    # colour handling
    cmap = linear_cmap("y", palette=palette, low=np.nanmin(y), high=np.nanmax(y))

    plot.scatter(x="time", y="y", source=source, color=cmap, marker="circle", legend_label=f"{lc.survey} {lc.band}")

    # plot errors
    err_xs = [[t, t] for t in time]
    err_ys = [[v - e, v + e] for v, e in zip(y, y_err)]
    err_source = ColumnDataSource(data=dict(xs=err_xs, ys=err_ys, y_val=y))
    err_cmap = linear_cmap("y_val", palette=palette, low=np.nanmin(y), high=np.nanmax(y))
    plot.multi_line(
        xs="xs", ys="ys", source=err_source, color=err_cmap, line_width=0.5, level="underlay", legend_label=f"{lc.survey} {lc.band}"
    )

    # don't show MJD in scientific notation
    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    return plot


def plot(lightcurves: list[Lightcurve], *args, **kwargs):
    """
    Plots any number of light curves into combined per-survey plots
    """

    colours = kwargs.get("colours")
    bands = kwargs.get("bands")
    time_format = kwargs.get("time_format", "reduced")

    plots = []
    surveys = list(set([lc.survey for lc in lightcurves]))
    print(surveys)

    # loop through surveys + combine light curve containers into single plot for each survey
    for survey in surveys:
        # check that all light curves share the same brightness type
        brightness_types = [lc.brightness_type for lc in lightcurves]
        if len(set(brightness_types)) > 1:
            raise ValueError("Invalid combination of lighcurve brightness types. Must be all 'flux' or all 'mag'.")

        # get single brightness type once above check has passed
        brightness_type = brightness_types[0]

        # filter data to only keep requested (and valid) light curves
        lcs = [lc for lc in lightcurves if lc.brightness_type and lc.survey == survey and (bands is None or lc.band in bands)]

        if not lcs:
            return None

        # colour handling
        palettes = assign_gradient_palettes(len(lcs), colours)

        # set up title
        band_names = ", ".join(d.band for d in lcs)

        # create per-survey plot
        plot = figure(
            width=400,
            height=400,
            title=f"{survey} {band_names} lightcurve(s)",
            x_axis_label="MJD" if time_format == "original" else "Time (days)",
            y_axis_label=brightness_type,
            tools=("pan,wheel_zoom,box_zoom,reset"),
        )

        # get MJD at start of data
        all_times = [t for lc in lcs for t in lc.mjd]
        time_min = min(all_times)

        # Plot each band independently
        for lc, palette in zip(lcs, palettes):
            plot = plot_band(plot=plot, lc=lc, palette=palette, time_min=time_min, time_format=time_format)

        if brightness_type == "flux":
            plot.y_range.flipped = True

        plots.append(plot)

    return [format_plot("lightcurve", p) for p in plots]
