import matplotlib
import numpy as np
from bokeh import events
from bokeh.models import BasicTickFormatter, ColumnDataSource, CustomJS
from bokeh.plotting import figure
from bokeh.transform import linear_cmap


class SupportedColours(object):
    def __init__(self, colour):
        self.colour = colour

    @property
    def get_cmap(self):
        if self.colour == "green":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["greenyellow", "forestgreen", "greenyellow"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "forestgreen"
        elif self.colour == "red":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["yellow", "red", "yellow"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "red"
        elif self.colour == "blue":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["aqua", "royalblue", "aqua"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "royalblue"
        elif self.colour == "black":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["lightgray", "black", "lightgray"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "black"
        elif self.colour == "orange":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["gold", "orange", "gold"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "orange"
        elif self.colour == "purple":
            colourmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "", ["orchid", "darkviolet", "orchid"]
            )
            palette = [
                matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0, 1, 255))
            ]
            error_colour = "darkviolet"
        else:
            raise ValueError(
                "Unsupported colour. Supported colours are: green, red, blue, black, orange, purple"
            )

        return palette, error_colour


def plot_data(plot, band, colour, time_min, survey, timeformat):
    if timeformat == "reduced":
        if using_original_times:
            time = [x - time_min for x in band[time_unit]]
        else:
            time = band[time_unit]
    elif timeformat == "original":
        time = band[time_unit]
    else:
        raise ValueError("Invalid timeformat. Accepted values: reduced, original.")

    palette, error_colour = SupportedColours(colour).get_cmap
    source = ColumnDataSource(data={"time": time, "mag": band["mag"]})
    cmap = linear_cmap(
        field_name="mag", palette=palette, low=min(band["mag"]), high=max(band["mag"])
    )

    plot.scatter(
        source=source,
        x="time",
        y="mag",
        color=cmap,
        marker="circle",
        legend_label=f"{survey} {band['band']}",
    )

    err_xs = [[x, x] for x in time]
    err_ys = [[y - y_err, y + y_err] for y, y_err in zip(band["mag"], band["mag_err"])]
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
    if bands:
        data = [x for x in struct.data if x["band"] in bands]
    else:
        data = [x for x in struct.data if x["mag"] is not None]

    if len(data) == 0:
        print("Note: Could not plot light curve, no data found.")
        return None

    if not colours:
        colours = ["black" for band in data]
    if len(colours) < len(data):
        for i in range(0, len(data) - len(colours)):
            colours.append("black")
    if not isinstance(colours, list):
        try:
            colours = [colours]
        except:
            raise Exception(
                f"Invalid colours type. Expected str or list, got {type(colours)}."
            )

    lightcurve_bands = ""
    for band in data:
        lightcurve_bands += f"{band['band']},"
    lightcurve_bands = lightcurve_bands.rstrip(",")

    plot = figure(
        width=400,
        height=400,
        title=f"{struct.survey} {lightcurve_bands} lightcurve(s)",
        x_axis_label="Observation Date / days",
        y_axis_label=f"{lightcurve_bands}",
    )

    global using_original_times, time_unit
    using_original_times = False

    if "hjd_ori" in band:
        time_unit = "hjd"
        if band["hjd_ori"]:
            time_unit += "_ori"
            using_original_times = True
    elif "mjd_ori" in band:
        time_unit = "mjd"
        if band["mjd_ori"]:
            time_unit += "_ori"
            using_original_times = True

    combined_times = []
    for band in data:
        combined_times += band[time_unit]

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
