import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, LinearColorMapper
from bokeh.models.formatters import BasicTickFormatter
from bokeh.plotting import figure

from ...structures.Lightcurve import Lightcurve
from ..colours import GRADIENT_MAPS, assign_gradient_palettes, get_gradient
from ..formatting import format_plot

MEAN_WARP_SCALE = 0.25


def plot_band(plot: figure, lc: Lightcurve, palette: list[str], time_min: float, time_format: str, cmap: str):
    """
    Plots a single light curves in a given band to an existing figure
    """

    # time handling
    time = lc.time
    if time_format == "reduced":
        time = [t - time_min for t in time]

    if lc.obj_id is not None:
        obj_id = str(lc.obj_id)
    else:
        obj_id = None

    # get brightness and error columns
    y = lc.brightness
    y_err = lc.brightness_err

    # calculate mean colour mapping
    mean_mag = np.nanmean(y)
    y_dev = np.abs(y - mean_mag)
    scale = np.nanmedian(y_dev) * MEAN_WARP_SCALE  # sets central band width
    y_warp = np.arcsinh(y_dev / scale)
    y_norm = y_warp / np.nanmax(y_warp)  # normalise

    # set up colour map
    colour_mapper = LinearColorMapper(palette=palette, low=0, high=1)

    df = pd.DataFrame({"time": time, "y": y, "y_err": y_err, "obj_id": obj_id})
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

    if lc.time_type == "mjd":
        legend_str = f"{lc.survey} {lc.band}"
    elif lc.time_type == "phase":
        legend_str = f"{lc.survey} {lc.band}\n{lc.fopt.value:.3f} {lc.fopt.unit.to_string('unicode')}"

    # plot errors
    err_xs = [[t, t] for t in time]
    err_ys = [[v - e, v + e] for v, e in zip(y, y_err)]
    err_source = ColumnDataSource(data={"xs": err_xs, "ys": err_ys, "y": y, "colour": y_norm})
    plot.multi_line(xs="xs", ys="ys", source=err_source, color=colour, line_width=0.5, legend_label=legend_str)

    # don't show MJD in scientific notation
    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    scatter = plot.scatter(x="time", y="y", source=source, color=colour, marker="circle", legend_label=legend_str)

    hvr.renderers = [scatter]
    if obj_id is not None:
        plot.add_tools(hvr)

    return plot


def dispatch_groups(lcs: list[Lightcurve], palette_map: dict, **kwargs: dict):
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

    # set up title
    band_names = ", ".join(d.band for d in lcs)

    # get x-axis label
    if getattr(lcs[0], "phase") is not None:
        x_label = "Phase"
    elif time_format == "original":
        x_label = "MJD"
    else:
        x_label = "Time (days)"

    # create plot(s)
    survey = lcs[0].survey

    multiband = list(set([lc.multiband for lc in lcs]))
    if len(multiband) > 1:
        raise ValueError("Plotting received mix of multiband and non-multiband data.")

    # get MJD at start of data
    if lcs[0].time_type == "mjd":
        all_times = [t for lc in lcs for t in lc.mjd]
        if not all_times:
            raise Exception("No time data found.")
        time_min = min(all_times)
    else:
        all_times = None
        time_min = None

    if multiband[0] is True or multiband[0] is None:
        plot = figure(
            width=400,
            height=400,
            title=f"{survey.upper()} {band_names} lightcurve(s)",
            x_axis_label=x_label,
            y_axis_label=brightness_type,
            tools=("pan,wheel_zoom,box_zoom,reset"),
        )

        # Plot each band independently
        for lc in lcs:
            plot = plot_band(plot=plot, lc=lc, palette=palette_map[lc.band], time_min=time_min, time_format=time_format, cmap=kwargs.get("cmap", "mean"))
        plots = [plot]

    elif multiband[0] is False:
        plots = []

        # Plot each band independently
        for lc in lcs:
            plot = figure(
                width=400,
                height=400,
                title=f"{survey.upper()} {band_names} lightcurve(s)",
                x_axis_label=x_label,
                y_axis_label=brightness_type,
                tools=("pan,wheel_zoom,box_zoom,reset"),
            )
            plot = plot_band(plot=plot, lc=lc, palette=palette_map[lc.band], time_min=time_min, time_format=time_format, cmap=kwargs.get("cmap", "mean"))
            plots.append(plot)

    if brightness_type != "flux":
        for plot in plots:
            plot.y_range.flipped = True

    return plots


def assign_band_colours(bands: list[str], colours: list[str] | None = None) -> dict[str, str]:
    cycle = [c for c in GRADIENT_MAPS if c != "black"]

    result = []

    if not colours:
        result = [cycle[i % len(cycle)] for i in range(len(bands))]
    else:
        colours = list(colours)

        # preserve user order FIRST
        result = colours[:]

        # fill remaining without duplicates
        used = list(dict.fromkeys(result))
        remaining = [c for c in cycle if c not in used]

        for c in remaining:
            if len(result) >= len(bands):
                break
            result.append(c)

        i = 0
        while len(result) < len(bands):
            result.append(cycle[i % len(cycle)])
            i += 1

    return {band: colour for band, colour in zip(bands, result)}


def get_band_colours(bands: list[str], colours: list[str] | None = None, gradient_size: int = 256, cycles: int = 1, reverse: bool = False) -> dict[str, list[str]]:
    band_colour_names = assign_band_colours(bands, colours)

    return {band: get_gradient(colour, gradient_size, cycles=cycles, reverse=reverse) for band, colour in band_colour_names.items()}


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


def plot(lightcurves: list[Lightcurve], *args: tuple, **kwargs: dict):
    """
    Plots any number of light curves into combined per-survey and per-id plots
    """

    plots = []

    # loop through surveys (will be removed)
    requested_bands = kwargs.get("bands")
    # filter data to only keep requested (and valid) light curves
    lcs = [lc for lc in lightcurves if lc.brightness_type and (requested_bands is None or lc.band in requested_bands)]
    # per_survey_lcs.sort(key=lambda lc: lc.obj_id)

    if kwargs.get("bands"):
        all_bands = kwargs["bands"]
        for band in kwargs["all_bands"]:
            if band in all_bands:
                continue
            all_bands.append(band)
    else:
        all_bands = kwargs["all_bands"]

    colours = kwargs.get("colours")
    palette_map = get_band_colours(all_bands, colours)

    # group by object ID
    lc_groups = group_lc_ids(lcs)
    for per_id_lcs in lc_groups:
        if len(set(lc.time_type for lc in per_id_lcs)) > 1:
            raise ValueError("Detected multiple time formats 'mjd' and 'phase' in Lightcurve plotting.")

        if per_id_lcs[0].time_type == "phase":
            kwargs["time_format"] = "original"
            kwargs["cmap"] = "flat"

        per_id_plots = dispatch_groups(per_id_lcs, palette_map, **kwargs)
        plots.extend(per_id_plots)

        if len(per_id_plots) == 1:
            for lc in per_id_lcs:
                lc._plot_id = per_id_plots[0].id
        else:
            for lc, plot in zip(per_id_lcs, per_id_plots):
                lc._plot_id = plot.id

    return [format_plot("lightcurve", p) for p in plots]
