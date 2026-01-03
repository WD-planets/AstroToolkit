import numpy as np
from bokeh.models import ColumnDataSource
from bokeh.models.formatters import BasicTickFormatter
from bokeh.plotting import figure
from bokeh.transform import linear_cmap

from ...structures.definitions import Lightcurve
from ..colours import GRADIENT_MAPS, get_gradient


def plot_band(plot: figure, lc: Lightcurve, colour, survey, brightness_type, time_min, time_format):
    # time handling
    time = lc.mjd
    if time_format == "reduced":
        time = [t - time_min for t in time]

    # get brightness and error columns
    y = getattr(lc, lc.brightness_type)
    y_err = getattr(lc, f"{brightness_type}_err")

    source = ColumnDataSource(data=dict(time=time, y=y, yerr=y_err))

    # colour handling
    palette, error_colour = get_gradient(colour)
    cmap = linear_cmap("y", palette=palette, low=np.nanmin(y), high=np.nanmax(y))

    plot.scatter(x="time", y="y", source=source, color=cmap, marker="circle", legend_label=f"{survey} {lc.band}")

    # plot errors
    err_xs = [[t, t] for t in time]
    err_ys = [[v - e, v + e] for v, e in zip(y, y_err)]
    plot.multi_line(err_xs, err_ys, color=error_colour, line_width=0.5, level="underlay")

    # don't show MJD in scientific notation
    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    return plot


def plot_lightcurve(lightcurves: list[Lightcurve], *args, **kwargs):
    """
    Plots any number of light curves into combined per-survey plots
    """

    colours = kwargs.get("colours")
    bands = kwargs.get("bands")
    time_format = kwargs.get("time_format", "reduced")

    plots = []
    surveys = [lc.survey for lc in lightcurves]

    # loop through surveys + combine light curve containers into single plot for each survey
    for survey in surveys:
        # check that all light curves share the same brightness type
        brightness_types = [lc.brightness_type for lc in lightcurves]
        if len(set(brightness_types)) > 1:
            raise ValueError("Invalid combination of lighcurve brightness types. Must be all 'flux' or all 'mag'.")

        # get single brightness type once above check has passed
        brightness_type = brightness_types[0]

        # filter data to only keep requested (and valid) light curves
        lightcurves = [lc for lc in lightcurves if lc.brightness_type and lc.survey == survey and (bands is None or lc.band in bands)]

        if not lightcurves:
            print("Note: No data to plot.")
            return None

        # colour handling
        available_colours = [c for c in GRADIENT_MAPS if c != "black"]
        if not colours:
            colours = ["black"] * len(lightcurves)
        elif len(colours) < len(lightcurves):
            fill = available_colours
            colours = colours + [fill[i % len(fill)] for i in range(len(lightcurves) - len(colours))]

        # set up title
        band_names = ", ".join(d["band"] for d in lightcurves)

        # create per-survey plot
        plot = figure(
            width=400,
            height=400,
            title=f"{survey} {band_names} lightcurve(s)",
            x_axis_label="MJD" if time_format == "original" else "Time (days)",
            y_axis_label=brightness_type,
        )

        # get MJD at start of data
        all_times = [t for lc in lightcurves for t in lc.mjd]
        time_min = min(all_times)

        # Plot each band independently
        for lc, colour in zip(lightcurves, colours):
            plot = plot_band(
                plot=plot, lc=lc, colour=colour, survey=survey, brightness_type=brightness_type, time_min=time_min, time_format=time_format
            )

        if brightness_type == "flux":
            plot.y_range.flipped = True

        plots.append(plot)

    return plots
