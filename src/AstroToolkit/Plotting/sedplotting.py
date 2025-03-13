import math

import numpy as np
from bokeh import events
from bokeh.models import (BasicTickFormatter, ColumnDataSource, CustomJS,
                          HoverTool)
from bokeh.plotting import figure


def plot_sed(struct, spectrum_overlay):
    colours = [
        "springgreen",
        "royalblue",
        "gold",
        "aquamarine",
        "deepskyblue",
        "orangered",
        "orange",
        "red",
        "black",
        "grey",
    ]

    plot = figure(
        width=400,
        height=400,
        title="Spectral Energy Distribution",
        x_axis_label="Effective Wavelength / \u212b",
        y_axis_label=r"\[\text{flux / mJy}\]",
        x_axis_type="log",
        y_axis_type="log",
    )

    plot.yaxis.formatter = BasicTickFormatter(use_scientific=False)
    plot.xaxis.major_label_overrides = {100000: r"\[10^5\]", 200000: r"\[2\times10^5\]"}

    plot.yaxis.ticker.desired_num_ticks = 5
    plot.xaxis.ticker.desired_num_ticks = 3

    if spectrum_overlay:
        survey = spectrum_overlay.survey
        spectrum_overlay = spectrum_overlay.data
        spectrum_x = spectrum_overlay["wavelength"]
        spectrum_y = spectrum_overlay["flux"]

        c = 2.988 * 10**18
        fnl = 1 * 10 ** (-23)
        spectrum_y_mjy = [y / ((fnl * c) / x**2) * 1000 for y, x in zip(spectrum_y, spectrum_x)]

        plot.line(spectrum_x, spectrum_y_mjy, color="black", line_width=1, legend_label=f"{survey} spectrum")

    legend_exists = False

    renderers = []
    for index, survey in enumerate(struct.data):
        for mag_name, wavelength, flux, flux_err in zip(
            survey["mag_name"], survey["wavelength"], survey["flux"], survey["flux_rel_err"]
        ):
            if math.isnan(flux_err):
                marker_type = "cross"
                suffix = "(Upper Limit)"
                size = 7.5
            else:
                marker_type = "circle"
                suffix = ""
                size = 5

            source = ColumnDataSource({"x": [wavelength], "y": [flux], "mag_name": [f"{survey['survey']} {mag_name}"]})

            hvr = HoverTool(tooltips=[("mag_name", "@mag_name"), ("wavelength", "@x")])
            scatter = plot.scatter(
                x="x",
                y="y",
                source=source,
                color=colours[index],
                marker=marker_type,
                legend_label=f"{survey['survey']} {suffix}",
                size=size,
            )

            renderers.append(scatter)

            err_xs = [[x, x] for x in [wavelength]]
            err_ys = [[y - y_err, y + y_err] for y, y_err in zip([flux], [flux_err])]
            if not np.isnan(flux_err):
                plot.multi_line(
                    err_xs,
                    err_ys,
                    color=colours[index],
                    legend_label=f"{survey['survey']} {suffix}",
                    line_width=0.5,
                    line_cap="square",
                )

            legend_exists = True

    hvr.renderers = renderers
    plot.add_tools(hvr)

    if legend_exists:
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

        plot.legend.location = "bottom_right"

    return plot
