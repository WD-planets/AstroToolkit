import math

from bokeh import events
from bokeh.models import BasicTickFormatter, CustomJS
from bokeh.plotting import figure


def plot_sed(struct, spectrum_overlay, survey):
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
        x_axis_label=r"\[\lambda_{\text{eff}}\text{ }[\text{AA}]\]",
        y_axis_label=r"\[\text{flux [mJy]}\]",
        x_axis_type="log",
        y_axis_type="log",
    )

    plot.yaxis.formatter = BasicTickFormatter(use_scientific=False)
    plot.xaxis.formatter = BasicTickFormatter(use_scientific=False)

    plot.xaxis.ticker.desired_num_ticks = 4

    if spectrum_overlay:
        if not survey:
            raise ValueError("Survey required for SED spectrum overlay.")

        from ..Tools import query

        if struct.source:
            spectrum_data = query(
                kind="spectrum", survey=survey, source=struct.source
            ).data
        else:
            spectrum_data = query(kind="spectrum", survey=survey, pos=struct.pos).data

        if spectrum_data:
            spectrum_x = spectrum_data["wavelength"]
            spectrum_y = spectrum_data["flux"]

            c = 2.988 * 10**18
            fnl = 1 * 10 ** (-23)
            spectrum_y_mjy = [
                y / ((fnl * c) / x**2) * 1000 for y, x in zip(spectrum_y, spectrum_x)
            ]

            plot.line(spectrum_x, spectrum_y_mjy, color="black", line_width=1)

    legend_exists = False
    for index, survey in enumerate(struct.data):
        for x, y, y_err in zip(
            survey["wavelength"], survey["flux"], survey["flux_rel_err"]
        ):
            if math.isnan(y_err):
                marker_type = "cross"
                suffix = "(Upper Limit)"
                size = 7.5
            else:
                marker_type = "circle"
                suffix = ""
                size = 5

            plot.scatter(
                x=x,
                y=y,
                color=colours[index],
                marker=marker_type,
                legend_label=f"{survey['survey']} {suffix}",
                size=size,
            )

            plot.multi_line(
                [x, x],
                [y - y_err, y + y_err],
                color=colours[index],
                legend_label=f"{survey['survey']} {suffix}",
                line_width=0.5,
                line_cap="square",
            )

            legend_exists = True

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

    return plot
