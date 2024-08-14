import numpy as np
from bokeh import events
from bokeh.models import CustomJS, Range1d
from bokeh.plotting import figure


def plot_image(struct):
    plot = figure(
        width=400, height=400, title=f'{struct.survey} Image ({struct.data["size"]}")'
    )

    plot.min_border = 75
    plot.grid.grid_line_color = None

    image_data = struct.data["image_data"]
    image_header = struct.data["image_header"]
    xlim, ylim = image_header["NAXIS1"], image_header["NAXIS2"]

    if xlim != ylim:
        xlim = ylim

    wcs = struct.data["wcs"]
    x_points, y_points = (
        np.arange(start=0, stop=xlim + 1, step=1),
        np.arange(start=0, stop=ylim + 1, step=1),
    )

    coords = wcs.all_pix2world(x_points, y_points, 1)
    x_points, y_points = coords[0], coords[1]

    x_range = max(x_points) - min(x_points)
    y_range = max(y_points) - min(y_points)

    plot.x_range = Range1d(max(x_points), min(x_points))
    plot.y_range = Range1d(min(y_points), max(y_points))

    plot.image(
        image=[image_data],
        x=x_points[0],
        y=y_points[0],
        dw=x_range,
        dh=y_range,
        palette="Greys256",
        level="image",
        origin="bottom_right",
        anchor="bottom_right",
    )

    focus_ra, focus_dec = struct.data["image_focus"]
    plot.scatter(x=focus_ra, y=focus_dec, marker="cross", color="black", size=10)

    if struct.data["overlay"]:
        pass
    else:
        return plot

    legend = False
    if "overlay" in struct.data:
        overlay_data = struct.data["overlay"]
        for data_point in overlay_data:
            if data_point["overlay_type"] == "detection_mag":
                if data_point["corrected"]:
                    line_style = "solid"
                    legend_suffix = "(corrected)"
                if not data_point["corrected"]:
                    line_style = "dotted"
                    legend_suffix = "(uncorrected)"
                if data_point["marker_type"] == "circle":
                    plot.circle(
                        x=data_point["ra"],
                        y=data_point["dec"],
                        radius=data_point["marker_size"],
                        line_color=data_point["colour"],
                        line_width=2,
                        line_dash=line_style,
                        fill_color=None,
                        legend_label=f"{data_point['survey']} {data_point['mag_name']} {legend_suffix}",
                    )
            elif data_point["overlay_type"] == "detection":
                if data_point["corrected"]:
                    marker = "cross"
                    legend_suffix = "(corrected)"
                if not data_point["corrected"]:
                    marker = "x"
                    legend_suffix = "(uncorrected)"
                plot.scatter(
                    x=data_point["ra"],
                    y=data_point["dec"],
                    marker=marker,
                    color=data_point["colour"],
                    legend_label=f"{data_point['survey']} {legend_suffix}",
                    size=20,
                )
            elif data_point["overlay_type"] == "tracer":
                plot.scatter(
                    x=data_point["ra"],
                    y=data_point["dec"],
                    marker="dot",
                    color=data_point["colour"],
                    size=20,
                    legend_label=f"{data_point['survey']} tracer",
                )
            legend = True

    if legend:
        plot.legend.click_policy = "hide"

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
