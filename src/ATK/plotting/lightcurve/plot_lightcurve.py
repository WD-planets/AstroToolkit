import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper
from bokeh.models.formatters import BasicTickFormatter
from bokeh.plotting import figure

from ...structures.definitions import Lightcurve
from ..colours import assign_gradient_palettes
from ..formatting import format_plot

MEAN_WARP_SCALE = 0.25


def plot_band(plot: figure, lc: Lightcurve, palette: list[str], time_min: float, time_format: str, cmap: str):
    """
    Plots a single light curves in a given band to an existing figure
    """

    # time handling
    time = lc.mjd
    if time_format == "reduced":
        time = [t - time_min for t in time]

    # get brightness and error columns
    obj_id = str(lc.obj_id)
    y = getattr(lc, lc.brightness_type)
    y_err = getattr(lc, f"{lc.brightness_type}_err")

    # calculate mean colour mapping
    mean_mag = np.nanmean(y)
    y_dev = np.abs(y - mean_mag)
    scale = np.nanmedian(y_dev) * MEAN_WARP_SCALE  # sets central band width
    y_warp = np.arcsinh(y_dev / scale)
    y_norm = y_warp / np.nanmax(y_warp)  # normalise

    # set up colour map
    colour_mapper = LinearColorMapper(palette=palette, low=0, high=1)

    df = pd.DataFrame(
        {"time": time, "y": getattr(lc, lc.brightness_type), "y_err": getattr(lc, f"{lc.brightness_type}_err"), "obj_id": obj_id}
    )
    df["colour"] = y_norm

    source = ColumnDataSource(df)
    hvr = HoverTool(tooltips=[("obj_id", "@obj_id")])

    #  final colours
    if cmap == "mean":
        colour = {"field": "colour", "transform": colour_mapper}
    elif cmap == "flat":
        colour = palette[0]
    else:
        raise ValueError(f"Invalid cmap '{cmap}'.")

    scatter = plot.scatter(x="time", y="y", source=source, color=colour, marker="circle", legend_label=f"{lc.survey} {lc.band}")

    hvr.renderers = [scatter]
    plot.add_tools(hvr)

    # plot errors
    err_xs = [[t, t] for t in time]
    err_ys = [[v - e, v + e] for v, e in zip(y, y_err)]
    err_source = ColumnDataSource(data={"xs": err_xs, "ys": err_ys, "y": y, "colour": y_norm})
    plot.multi_line(
        xs="xs", ys="ys", source=err_source, color=colour, line_width=0.5, level="underlay", legend_label=f"{lc.survey} {lc.band}"
    )

    # don't show MJD in scientific notation
    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    return plot


def group_lc_ids(lcs: list[Lightcurve]):
    """
    Group light curves by object ID
    """

    ids = list(set([lc.obj_id for lc in lcs]))

    groups = []
    for id in ids:
        id_group = [lc for lc in lcs if lc.obj_id == id]
        groups.append(id_group)

    return groups


def dispatch_groups(survey: str, lcs: list[Lightcurve], **kwargs: dict):
    """
    Plots light curves of a given survey grouped by object ID
    """

    time_format = kwargs.get("time_format", "reduced")

    # check that all light curves share the same brightness type
    brightness_types = [lc.brightness_type for lc in lcs]
    if len(set(brightness_types)) > 1:
        raise ValueError("Invalid combination of lighcurve brightness types. Must be all 'flux' or all 'mag'.")

    # get single brightness type once above check has passed
    brightness_type = brightness_types[0]

    # colour handling
    palettes = assign_gradient_palettes(len(lcs), kwargs.get("colours"))

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
        plot = plot_band(plot=plot, lc=lc, palette=palette, time_min=time_min, time_format=time_format, cmap=kwargs.get("cmap", "mean"))

    if brightness_type == "flux":
        plot.y_range.flipped = True

    return plot


def plot(lightcurves: list[Lightcurve], *args: tuple, **kwargs: dict):
    """
    Plots any number of light curves into combined per-survey and per-id plots
    """

    bands = kwargs.get("bands")

    plots = []
    surveys = list(set([lc.survey for lc in lightcurves]))

    # loop through surveys + combine light curve containers into single plot for each survey
    for survey in surveys:
        # filter data to only keep requested (and valid) light curves
        per_survey_lcs = [lc for lc in lightcurves if lc.brightness_type and lc.survey == survey and (bands is None or lc.band in bands)]

        if not per_survey_lcs:
            continue

        # group by object ID
        lc_groups = group_lc_ids(per_survey_lcs)
        for per_id_lcs in lc_groups:
            per_id_plots = dispatch_groups(survey, per_id_lcs, **kwargs)
            plots.append(per_id_plots)

    return [format_plot("lightcurve", p) for p in plots]
