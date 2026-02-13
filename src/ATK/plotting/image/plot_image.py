import astropy.units as u
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.wcs.utils import proj_plane_pixel_scales
from bokeh.models import (ColumnDataSource, HoverTool, LinearColorMapper,
                          OpenURL, Range1d, TapTool)
from bokeh.palettes import Greys256, Viridis256
from bokeh.plotting import figure

from ...configuration.base_config import BASE_CONFIG
from ...structures.Image import Image
from ...utilities.coordinates import correct_dataframe_coords
from ..colours import get_palette
from ..formatting import format_plot
from .false_colour import get_false_cmap


def get_simbad_urls(image: Image, overlay_data: pd.DataFrame) -> pd.DataFrame:
    """
    Adds simbad_url (+ simbad_ra/simbad_dec) columns to overlay dataframe
    """

    # get j2000 coords of detections + add to "simbad_ra" / "simbad_dec" columns
    overlay_data = correct_dataframe_coords(
        overlay_data, image.epoch, Time("2000-01-01", format="iso"), output_cols=["simbad_ra", "simbad_dec"]
    )

    simbad_radius = BASE_CONFIG.get("overlay_settings", "simbad_radius")

    # add url column
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


def get_marker_size(image: Image, overlay_data: pd.DataFrame, relative_axes: bool) -> float:
    """
    Adds marker_size (for scaling detection circles) and pointer_size (central detection pointer) columns to the overlay dataframe
    """

    # basic constants determined by-eye
    BASE_RADIUS = 5e-4
    SCALE_FACTOR = 5e-3
    HALF_IMAGE = image.size / 3600
    POINTER_SIZE_RATIO = 1 / 7.5

    # somewhat physically-based scaling relation (compare to magnitude 20, logarithmic scaling with flux)
    flux_multiplier = 10 ** (-0.4 * (overlay_data["mag"] - 20))

    # np.log10() dampens the effect on brighter sources so that marker sizes don't explode
    marker_sizes = BASE_RADIUS + np.log10(1 + flux_multiplier) * HALF_IMAGE * SCALE_FACTOR

    # conversion to arcsecond-based marker sizes
    if relative_axes:
        marker_sizes *= 3600

    overlay_data["marker_radius"] = marker_sizes
    overlay_data["pointer_radius"] = marker_sizes * POINTER_SIZE_RATIO

    return overlay_data


def plot_overlay(plot: figure, image: Image, relative_axes: bool) -> figure:
    """
    Adds an overlay of nearby detections to the image
    """

    # make copy so that original overlay dataframe isn't edited
    overlay = image.overlay.copy()

    # set up hovertool
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

    # set up taptool
    taptool = TapTool(renderers=hvr.renderers, callback=OpenURL(url="@simbad_url"))
    taptool.renderers = []

    # set up label column for legend entries
    overlay["label"] = overlay["survey"].astype(str)
    overlay.loc[overlay["mag"].isna(), "label"] += " detection"
    overlay.loc[overlay["mag"].notna(), "label"] += " " + overlay.loc[overlay["mag"].notna(), "mag_name"]
    overlay["gaia_match"] = overlay["gaia_match"].astype(str)

    # get simbad URLs for taptool
    overlay = get_simbad_urls(image, overlay)

    # set marker locations
    if relative_axes:
        overlay["marker_ra"] = (overlay["ra"] - image.search_pos.ra.value) * 3600
        overlay["marker_dec"] = (overlay["dec"] - image.search_pos.dec.value) * 3600
    else:
        overlay["marker_ra"] = overlay["ra"]
        overlay["marker_dec"] = overlay["dec"]

    magnitudes = sorted(overlay["mag_name"].unique())
    # get unique colours for markers, shift = 1 avoids blue markers for viridis
    colours = get_palette(len(magnitudes), shift=1)
    # get dict of mag_name: colour and map to overlay dataframe
    colour_map = dict(zip(magnitudes, colours))
    overlay["colour"] = overlay["mag_name"].map(colour_map)

    # split into detections and magnitude-scaled detections
    mask = overlay["mag"].isna()
    non_nan_mag = overlay.loc[~mask].copy()
    nan_mag = overlay.loc[mask].copy()

    # get marker sizes
    non_nan_mag = get_marker_size(image, non_nan_mag, relative_axes)

    # magnitude-scaled detections
    for (survey, label), group in non_nan_mag.groupby(["survey", "label"]):
        plot.circle(
            source=ColumnDataSource(group),
            x="marker_ra",
            y="marker_dec",
            radius="marker_radius",
            line_color="colour",
            line_width=2,
            fill_color=None,
            legend_label=label,
        )

        clickable_marker = plot.circle(
            source=ColumnDataSource(group),
            x="marker_ra",
            y="marker_dec",
            radius="pointer_radius",
            line_color="colour",
            line_width=2,
            fill_color="colour",
            alpha=0.5,
            legend_label=label,
        )
        hvr.renderers.append(clickable_marker)
        taptool.renderers.append(clickable_marker)

    # non-scaled detections
    for (survey, label), group in nan_mag.groupby(["survey", "label"]):
        scatter = plot.scatter(
            source=ColumnDataSource(group),
            x="marker_ra",
            y="marker_dec",
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
    """
    Plots an ATK Image object
    """

    # create figure
    plot = figure(width=400, height=400, title=f"{image.survey} {image.band}-band Image ({image.size})", tools=("pan,wheel_zoom,reset"))
    plot.grid.grid_line_color = None

    # get image centre and deg/pixel
    n_pixels = (image.hdu.data.shape[1], image.hdu.data.shape[0])
    pixel_scales = proj_plane_pixel_scales(image.wcs)

    # get raw data array
    image_data = image.hdu.data

    # horizontally flip image for Bokeh plotting
    image_data = np.fliplr(image_data)

    # subtract min, asin stretch + percentile squash image
    nan_mask = np.isnan(image_data)
    image_data[nan_mask] = 0.0
    image_data = np.arcsinh(image_data - np.min(image_data))
    vmin, vmax = np.percentile(image_data, [5.0, 99.9])
    image_data[nan_mask] = np.nan

    # get colour map
    cmap = kwargs.get("cmap", "viridis")
    if cmap == "viridis":
        colour_mapper = LinearColorMapper(palette=Viridis256, low=vmin, high=vmax)
    elif cmap == "grey":
        colour_mapper = LinearColorMapper(palette=Greys256, low=vmin, high=vmax)
    elif cmap == "false_colour":
        colour_mapper = LinearColorMapper(palette=get_false_cmap(image.survey, image.band), low=vmin, high=vmax)

    relative_axes = kwargs.get("relative_axes", True)
    image_size_as = image.size.to(u.arcsec)

    # relative (+- arcsec from centre) axes
    if relative_axes:
        plot.xaxis.axis_label = "Relative Right Ascension / arcsec"
        plot.yaxis.axis_label = "Relative Declination / arcsec"

        # image bounds in arcseconds
        x_bounds = (-n_pixels[0] / 2 * pixel_scales[0] * 3600, n_pixels[0] / 2 * pixel_scales[0] * 3600)
        y_bounds = (-n_pixels[1] / 2 * pixel_scales[1] * 3600, n_pixels[1] / 2 * pixel_scales[1] * 3600)

        # image size in arcseconds
        x_range = (x_bounds[1] - x_bounds[0]) * u.arcsec
        y_range = (y_bounds[1] - y_bounds[0]) * u.arcsec

        # cap axes to image extent
        if x_range < image_size_as:
            plot.x_range = Range1d(x_bounds[1], x_bounds[0])
        else:
            plot.x_range = Range1d(image_size_as.value / 2, -image_size_as.value / 2)

        if y_range < image.size:
            plot.y_range = Range1d(y_bounds[0], y_bounds[1])
        else:
            plot.y_range = Range1d(-image_size_as.value / 2, image_size_as.value / 2)

        search_pos_ra, search_pos_dec = 0.0, 0.0

    # coordinates axes
    else:
        image_search_pos = (image.search_pos.ra.value, image.search_pos.dec.value)
        image_size_deg = image_size_as.to(u.deg)

        plot.xaxis.axis_label = "Right Ascension / deg"
        plot.yaxis.axis_label = "Declination / deg"

        # image bounds in deg
        x_bounds = (image_search_pos[0] - n_pixels[0] / 2 * pixel_scales[0], image_search_pos[0] + n_pixels[0] / 2 * pixel_scales[0])
        y_bounds = (image_search_pos[1] - n_pixels[1] / 2 * pixel_scales[1], image_search_pos[1] + n_pixels[1] / 2 * pixel_scales[1])

        # image range in deg
        x_range = (x_bounds[1] - x_bounds[0]) * u.deg
        y_range = (y_bounds[1] - y_bounds[0]) * u.deg

        # cap axes to image extent
        if x_range < image_size_deg:
            plot.x_range = Range1d(x_bounds[1], x_bounds[0])
        else:
            plot.x_range = Range1d(image_search_pos[0] + image_size_deg.value / 2, image_search_pos[0] - image_size_deg.value / 2)

        if y_range < image_size_deg:
            plot.y_range = Range1d(y_bounds[0], y_bounds[1])
        else:
            plot.y_range = Range1d(image_search_pos[1] - image_size_deg.value / 2, image_search_pos[1] + image_size_deg.value / 2)

        search_pos_ra, search_pos_dec = image_search_pos[0], image_search_pos[1]

    # don't allow panning outside of image bounds
    plot.x_range.bounds = "auto"
    plot.y_range.bounds = "auto"

    # plot image
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

    # plot search_pos marker
    plot.scatter(x=search_pos_ra, y=search_pos_dec, marker="cross", color="lime", size=25, line_width=4)

    # plot overlay
    if image.overlay is not None:
        plot = plot_overlay(plot, image, relative_axes)

    image._plot_id = plot.id

    return format_plot("image", plot)
