import matplotlib
import numpy as np
from bokeh import events
from bokeh.models import BasicTickFormatter, ColumnDataSource, CustomJS
from bokeh.plotting import figure
from bokeh.transform import linear_cmap

from ..Utility import getBrightnessType


class SupportedColours(object):
    def __init__(self, colour):
        self.colour = colour
        self.supported_colours = ["green", "red", "blue", "black", "orange", "purple"]

    @property
    def get_cmap(self):
        if self.colour == "green":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["greenyellow", "forestgreen", "greenyellow"]
            )
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "forestgreen"
        elif self.colour == "red":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["yellow", "red", "yellow"])
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "red"
        elif self.colour == "blue":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["aqua", "royalblue", "aqua"])
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "royalblue"
        elif self.colour == "black":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgray", "black", "lightgray"])
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "black"
        elif self.colour == "orange":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["gold", "orange", "gold"])
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "orange"
        elif self.colour == "purple":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["orchid", "darkviolet", "orchid"])
            palette = [matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))]
            error_colour = "darkviolet"

        return palette, error_colour


def plot_data(plot, band, colour, time_min, survey, timeformat):
    if timeformat == "reduced":
        time = [x - time_min for x in band["mjd"]]
    else:
        time = band["mjd"]

    palette, error_colour = SupportedColours(colour).get_cmap
    source = ColumnDataSource(data={"time": time, brightness_type: band[brightness_type]})
    cmap = linear_cmap(
        field_name=brightness_type, palette=palette, low=min(band[brightness_type]), high=max(band[brightness_type])
    )

    plot.scatter(
        source=source, x="time", y=brightness_type, color=cmap, marker="circle", legend_label=f"{survey} {band['band']}"
    )

    err_xs = [[x, x] for x in time]
    err_ys = [[y - y_err, y + y_err] for y, y_err in zip(band[brightness_type], band[f"{brightness_type}_err"])]
    plot.multi_line(
        err_xs,
        err_ys,
        color=error_colour,
        legend_label=f"{survey} {band['band']}",
        level="underlay",
        line_width=0.5,
        line_cap="square",
    )

    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    return plot


def plot_lightcurve(struct, colours, bands, timeformat):
    global brightness_type
    brightness_type = getBrightnessType(struct.data)

    if bands:
        data = [x for x in struct.data if x["band"] in bands and x[brightness_type] is not None]
    else:
        data = [x for x in struct.data if x[brightness_type] is not None]

    if len(data) == 0:
        print("Note: Could not plot light curve, no data found in requested bands.")
        return None

    if not colours:
        colours = ["black" for band in data]
    if len(colours) < len(data):
        for i in range(0, len(data) - len(colours)):
            colours.append("black")

    lightcurve_bands = ""
    for band in data:
        lightcurve_bands += f"{band['band']}, "
    lightcurve_bands = lightcurve_bands.rstrip(", ")
    if len(data) > 1:
        if brightness_type == "flux":
            lightcurve_bands += " fluxes"
        else:
            lightcurve_bands += " mags"
    elif len(data) == 1:
        if brightness_type == "flux":
            lightcurve_bands += " flux"
        else:
            lightcurve_bands += " mag"

    plot = figure(
        width=400,
        height=400,
        title=f"{struct.survey} {lightcurve_bands} lightcurve(s)",
        x_axis_label="Observation Date / days",
        y_axis_label=f"{lightcurve_bands}",
    )

    if timeformat == "original":
        plot.xaxis.axis_label = "MJD"

    combined_times = []
    for band in data:
        combined_times += band["mjd"]

    time_min = min(combined_times)

    for band, colour in zip(data, colours):
        plot = plot_data(plot, band, colour, time_min, struct.survey, timeformat)

    plot.y_range.flipped = True
    plot.legend.click_policy = "hide"

    # Double click to hide legend
    toggle_legend_js = CustomJS(
        args=dict(leg=plot.legend[0]),
        code="""
            if (leg.visible) {
                leg.visible = false
                }
            else {
                leg.visible = true
            }
    """,
    )

    plot.js_on_event(events.DoubleTap, toggle_legend_js)

    return plot
