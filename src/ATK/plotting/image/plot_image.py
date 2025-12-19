import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.wcs.utils import proj_plane_pixel_scales
from bokeh.models import (ColumnDataSource, HoverTool, LinearColorMapper,
                          OpenURL, Range1d, TapTool)
from bokeh.palettes import Greys256, Viridis256
from bokeh.plotting import figure

from ...configuration.base_config import BASE_CONFIG
from ...structures.definitions import Image
from ...utilities.coordinates import correct_dataframe_coords
from ..colours import get_palette
from ..formatting import format_plot
from .false_colour import get_false_cmap


def get_relative_axes(image: Image):
    wcs = image.wcs
    image_data = image.hdu.data
    n_pix_y, n_pix_x = image_data.shape

    # center pixel
    x_centre = (n_pix_x) / 2
    y_centre = (n_pix_y) / 2

    # pixel offsets from center
    x_pix = np.arange(n_pix_x + 1) - x_centre
    y_pix = np.arange(n_pix_y + 1) - y_centre

    # pixel scales
    pixscale_x, pixscale_y = proj_plane_pixel_scales(wcs)

    # Convert to arcsec
    x_arcsec_edges = x_pix * pixscale_x * 3600.0
    y_arcsec_edges = y_pix * pixscale_y * 3600.0

    return x_arcsec_edges, y_arcsec_edges


def get_simbad_urls(image: Image, overlay_data: pd.DataFrame):
    overlay_data = correct_dataframe_coords(
        overlay_data, image.epoch, Time("2000-01-01", format="iso"), output_cols=["simbad_ra", "simbad_dec"]
    )

    simbad_radius = BASE_CONFIG.get("overlay_settings", "simbad_radius")

    overlay_data["simbad_url"] = (
        "https://simbad.cds.unistra.fr/simbad/sim-coo?Coord="
        + overlay_data["simbad_ra"].astype(str)
        + "+"
        + overlay_data["simbad_dec"].astype(str)
        + "&CooFrame=icrs&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius="
        + str(simbad_radius)
        + "&Radius.unit=arcsec&submit=submit+query"
    )

    return overlay_data


def get_marker_size(image: Image, overlay_data: pd.DataFrame, relative_axes: bool):
    BASE_RADIUS = 5e-4
    SCALE_FACTOR = 5e-3
    HALF_IMAGE = image.size / 3600

    # somewhat physically-based scaling relation (compare to magnitude 20, logarithmic scaling with flux)
    flux_multiplier = 10 ** (-0.4 * (overlay_data["mag"] - 20))

    # np.log10() dampens the effect on brighter sources so that marker sizes don't explode
    marker_sizes = BASE_RADIUS + np.log10(1 + flux_multiplier) * HALF_IMAGE * SCALE_FACTOR

    if relative_axes:
        marker_sizes *= 3600

    overlay_data["marker_radius"] = marker_sizes
    overlay_data["pointer_radius"] = marker_sizes / 7.5

    return overlay_data


def plot_overlay(plot: figure, image: Image, relative_axes: bool):
    overlay = image.overlay

    hvr = HoverTool(
        tooltips=[
            ("survey", "@survey"),
            ("ra", "@ra"),
            ("dec", "@dec"),
            ("band", "@mag_name"),
            ("mag", "@mag"),
            ("error", "@err"),
            ("simbad_id", "@simbad_id"),
            ("corrected", "@gaia_match"),
        ]
    )
    hvr.renderers = []

    taptool = TapTool(renderers=hvr.renderers, callback=OpenURL(url="@simbad_url"))
    taptool.renderers = []

    overlay["label"] = overlay["survey"].astype(str)
    overlay.loc[overlay["mag"].isna(), "label"] += " detection"
    overlay.loc[overlay["mag"].notna(), "label"] += " " + overlay.loc[overlay["mag"].notna(), "mag_name"]
    overlay["gaia_match"] = overlay["gaia_match"].astype(str)

    overlay = get_simbad_urls(image, overlay)

    if relative_axes:
        overlay["ra"] = (overlay["ra"] - image.focus.ra.value) * 3600
        overlay["dec"] = (overlay["dec"] - image.focus.dec.value) * 3600

    magnitudes = sorted(overlay["mag_name"].unique())
    colours = get_palette(len(magnitudes), shift=1)
    colour_map = dict(zip(magnitudes, colours))
    overlay["colour"] = overlay["mag_name"].map(colour_map)

    mask = overlay["mag"].isna()
    non_nan_mag = overlay.loc[~mask].copy()
    nan_mag = overlay.loc[mask].copy()

    non_nan_mag = get_marker_size(image, non_nan_mag, relative_axes)

    for (survey, label), group in non_nan_mag.groupby(["survey", "label"]):
        plot.circle(
            source=ColumnDataSource(group),
            x="ra",
            y="dec",
            radius="marker_radius",
            line_color="colour",
            line_width=2,
            fill_color=None,
            legend_label=label,
        )

        clickable_marker = plot.circle(
            source=ColumnDataSource(group),
            x="ra",
            y="dec",
            radius="pointer_radius",
            line_color="colour",
            line_width=2,
            fill_color="colour",
            alpha=0.5,
            legend_label=label,
        )
        hvr.renderers.append(clickable_marker)
        taptool.renderers.append(clickable_marker)

    for (survey, label), group in nan_mag.groupby(["survey", "label"]):
        scatter = plot.scatter(
            source=ColumnDataSource(group),
            x="ra",
            y="dec",
            size=20,
            line_width=4,
            color="colour",
            marker="x",
            legend_label=label,
        )
        hvr.renderers.append(scatter)
        taptool.renderers.append(scatter)

    plot.add_tools(hvr)
    plot.add_tools(taptool)

    return plot


def plot(image: Image, *args: any, **kwargs: any) -> figure:
    plot = figure(
        width=400,
        height=400,
        title=f'{image.survey} {image.band}-band Image ({image.size}")',
        tools=("pan,wheel_zoom,reset"),
    )
    plot.grid.grid_line_color = None

    n_pixels = (image.hdu.data.shape[1], image.hdu.data.shape[0])
    image_focus = (image.focus.ra.value, image.focus.dec.value)
    pixel_scales = proj_plane_pixel_scales(image.wcs)

    # get raw data array
    image_data = image.hdu.data

    # horizontally flip image for Bokeh plotting
    image_data = np.fliplr(image_data)

    # subtract min, asin stretch + percentile squash
    nan_mask = np.isnan(image_data)
    image_data[nan_mask] = 0.0
    image_data = np.arcsinh(image_data - np.min(image_data))
    vmin, vmax = np.percentile(image_data, [5.0, 99.9])
    image_data[nan_mask] = np.nan

    cmap = kwargs.get("cmap", "viridis")
    if cmap == "viridis":
        colour_mapper = LinearColorMapper(palette=Viridis256, low=vmin, high=vmax)
    elif cmap == "grey":
        colour_mapper = LinearColorMapper(palette=Greys256, low=vmin, high=vmax)
    elif cmap == "false_colour":
        colour_mapper = LinearColorMapper(palette=get_false_cmap(image.survey, image.band), low=vmin, high=vmax)

    relative_axes = kwargs.get("relative_axes", True)

    if relative_axes:
        plot.xaxis.axis_label = "Relative Right Ascension / arcsec"
        plot.yaxis.axis_label = "Relative Declination / arcsec"

        x_bounds = (-n_pixels[0] / 2 * pixel_scales[0] * 3600, n_pixels[0] / 2 * pixel_scales[0] * 3600)
        y_bounds = (-n_pixels[1] / 2 * pixel_scales[1] * 3600, n_pixels[1] / 2 * pixel_scales[1] * 3600)

        x_range = x_bounds[1] - x_bounds[0]
        y_range = y_bounds[1] - y_bounds[0]

        if x_range < image.size:
            plot.x_range = Range1d(x_bounds[1], x_bounds[0])
        else:
            plot.x_range = Range1d(image.size / 2, -image.size / 2)

        if y_range < image.size:
            plot.y_range = Range1d(y_bounds[0], y_bounds[1])
        else:
            plot.y_range = Range1d(-image.size / 2, image.size / 2)

        focus_ra, focus_dec = 0.0, 0.0
    else:
        plot.xaxis.axis_label = "Right Ascension / deg"
        plot.yaxis.axis_label = "Declination / deg"

        x_bounds = (
            image_focus[0] - n_pixels[0] / 2 * pixel_scales[0],
            image_focus[0] + n_pixels[0] / 2 * pixel_scales[0],
        )
        y_bounds = (
            image_focus[1] - n_pixels[1] / 2 * pixel_scales[1],
            image_focus[1] + n_pixels[1] / 2 * pixel_scales[1],
        )

        x_range = x_bounds[1] - x_bounds[0]
        y_range = y_bounds[1] - y_bounds[0]

        if x_range < image.size:
            plot.x_range = Range1d(x_bounds[1], x_bounds[0])
        else:
            plot.x_range = Range1d(image_focus[0] + image.size / 2, image_focus[0] - image.size / 2)

        if y_range < image.size:
            plot.y_range = Range1d(y_bounds[0], y_bounds[1])
        else:
            plot.x_range = Range1d(image_focus[1] - image.size / 2, image_focus[1] + image.size / 2)

        focus_ra, focus_dec = image_focus[0], image_focus[1]

    plot.x_range.bounds = "auto"
    plot.y_range.bounds = "auto"

    plot.image(
        image=[image_data],
        x=x_bounds[0],
        y=y_bounds[0],
        dw=x_bounds[1] - x_bounds[0],
        dh=y_bounds[1] - y_bounds[0],
        level="image",
        origin="bottom_left",
        anchor="bottom_left",
        color_mapper=colour_mapper,
    )

    plot.scatter(x=focus_ra, y=focus_dec, marker="cross", color="lime", size=25, line_width=4)

    if image.overlay is not None:
        plot = plot_overlay(plot, image, relative_axes)

    return format_plot("image", plot)
