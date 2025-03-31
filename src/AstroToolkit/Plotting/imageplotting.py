import matplotlib.pyplot as plt
import numpy as np
from bokeh import events
from bokeh.models import (ColumnDataSource, CustomJS, HoverTool,
                          LinearColorMapper, NumeralTickFormatter, OpenURL,
                          Range1d, TapTool)
from bokeh.plotting import figure

from ..Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()


CMAP_MIN_VALUE = 25
CMAP_MAX_VALUE = 240


def plot_image(struct, simbad_search_radius=None):
    plot = figure(
        width=400,
        height=400,
        title=f'{struct.survey} Image ({struct.data["size"]}")',
        x_axis_label="Right Ascension / deg",
        y_axis_label="Declination / deg",
        tools=("pan,wheel_zoom,reset,tap"),
    )

    plot.grid.grid_line_color = None

    image_data = struct.data["image_data"]
    image_header = struct.data["image_header"]
    xlim, ylim = image_header["NAXIS1"], image_header["NAXIS2"]

    if xlim != ylim:
        xlim = ylim

    wcs = struct.data["wcs"]
    x_points, y_points = (np.arange(start=0, stop=xlim + 1, step=1), np.arange(start=0, stop=ylim + 1, step=1))

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
    plot.scatter(x=focus_ra, y=focus_dec, marker="cross", color="black", size=10, line_width=2)

    if struct.data["overlay"]:
        pass
    else:
        return plot

    legend = False

    clickable_markers = []
    if "overlay" in struct.data:
        overlay_data = struct.data["overlay"]
        colour_indices = [x["colour_index"] for x in overlay_data]
        index_range = max(colour_indices) - min(colour_indices)
        if index_range > 0:
            cmap_increment = (CMAP_MAX_VALUE - CMAP_MIN_VALUE) / index_range
        else:
            cmap_increment = 0
        for data_point in overlay_data:
            cmap = LinearColorMapper(palette="Turbo256", low=0, high=255)
            mapped_index = (data_point["colour_index"] - min(colour_indices)) * cmap_increment + CMAP_MIN_VALUE
            colour = {"field": "colour_index", "transform": cmap}
            source = ColumnDataSource(
                {"ra": [data_point["ra"]], "dec": [data_point["dec"]], "colour_index": [mapped_index]}
            )
            if data_point["overlay_type"] == "detection_mag":
                legend_label = f"{data_point['survey']} {data_point['mag_name']}"
                if data_point["marker_type"] == "circle":
                    plot.circle(
                        source=source,
                        x="ra",
                        y="dec",
                        radius=data_point["marker_size"],
                        line_color=colour,
                        line_width=3,
                        fill_color=None,
                        legend_label=legend_label,
                    )
            elif data_point["overlay_type"] == "detection":
                legend_label = f"{data_point['survey']}"
                plot.scatter(
                    source=source,
                    x="ra",
                    y="dec",
                    marker="cross",
                    color=colour,
                    legend_label=legend_label,
                    size=20,
                    line_width=3,
                )
            elif data_point["overlay_type"] == "tracer":
                legend_label = f"{data_point['survey']} tracer"
                plot.scatter(
                    x=data_point["ra"],
                    y=data_point["dec"],
                    marker="dot",
                    color=colour,
                    size=30,
                    legend_label=legend_label,
                )

            if data_point["overlay_type"] in ["detection_mag", "detection"]:
                from ..Tools import correctpm

                image_time = struct.data["image_time"]

                ra, dec = data_point["ra"], data_point["dec"]
                if data_point["corrected"]:
                    url_ra, url_dec = correctpm(
                        input_time=image_time,
                        target_time=[2000, 0],
                        pos=[ra, dec],
                        pmra=data_point["correction_pmra"],
                        pmdec=data_point["correction_pmdec"],
                    )
                else:
                    url_ra, url_dec = ra, dec

                url = f"https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={url_ra}+{url_dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={simbad_search_radius}&Radius.unit=arcsec&submit=submit+query"

                clickable_markers_source = ColumnDataSource(
                    {
                        "survey": [data_point["survey"]],
                        "ra": [data_point["ra"]],
                        "dec": [data_point["dec"]],
                        "mag": [data_point["mag"]],
                        "url_ra": [url_ra],
                        "url_dec": [url_dec],
                        "obj_id": [data_point["obj_id"]],
                        "simbad_id": [data_point["simbad_id"]],
                        "corrected": [data_point["corrected"]],
                        "url": [url],
                        "colour_index": [mapped_index],
                    }
                )

                clickable_marker = plot.circle(
                    source=clickable_markers_source,
                    x="ra",
                    y="dec",
                    radius=data_point["marker_size"] / 7.5,
                    line_color=colour,
                    line_width=2,
                    fill_color=colour,
                    alpha=0.5,
                    legend_label=legend_label,
                )
                clickable_markers.append(clickable_marker)

                legend = True

    hvr = HoverTool(
        tooltips=[
            ("survey", "@survey"),
            ("ra", "@ra"),
            ("dec", "@dec"),
            ("mag", "@mag"),
            ("survey_id", "@obj_id"),
            ("simbad_id", "@simbad_id"),
            ("corrected", "@corrected"),
        ]
    )
    hvr.renderers = clickable_markers
    plot.add_tools(hvr)

    taptool = plot.select(type=TapTool)[0]
    taptool.renderers = clickable_markers
    taptool.callback = OpenURL(url="@url")
    plot.add_tools(taptool)

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

    plot.xaxis.formatter = NumeralTickFormatter(format="0.000")
    plot.xaxis.ticker.desired_num_ticks = 3
    plot.yaxis.formatter = NumeralTickFormatter(format="0.000")
    plot.yaxis.ticker.desired_num_ticks = 3

    return plot
