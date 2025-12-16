import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales
from bokeh.models import LinearColorMapper, Range1d
from bokeh.palettes import Greys256, Viridis256
from bokeh.plotting import figure

from ...structures.definitions import Image
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


def plot(image: Image, *args: any, **kwargs: any) -> figure:
    plot = figure(
        width=400,
        height=400,
        title=f'{image.survey} {image.band}-band Image ({image.size}")',
        x_axis_label="Right Ascension / deg",
        y_axis_label="Declination / deg",
        tools=("pan,wheel_zoom,reset,tap"),
    )
    plot.grid.grid_line_color = None

    n_pixels = (image.hdu.data.shape[1], image.hdu.data.shape[0])
    image_focus = (image.focus.ra.value, image.focus.dec.value)
    pixel_scales = proj_plane_pixel_scales(image.wcs)

    # get raw data array
    image_data = image.hdu.data

    # horizontally flip image
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

    if kwargs.get("relative_axes", True):
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

    return format_plot("image", plot)
