import numpy as np
from astropy.wcs.utils import proj_plane_pixel_scales
from bokeh.models import LinearColorMapper, Range1d
from bokeh.palettes import Greys256, Viridis256
from bokeh.plotting import figure

from ...structures.definitions import Image
from ...utilities.defaults import EFFECTIVE_WAVELENGTHS
from ..formatting import format_plot
from .false_colour import wavelength_to_cmap


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

    xlim, ylim = image.hdu.header["NAXIS1"], image.hdu.header["NAXIS2"]

    if xlim != ylim:
        xlim = ylim

    x_points, y_points = (np.arange(start=0, stop=xlim + 1, step=1), np.arange(start=0, stop=ylim + 1, step=1))

    coords = image.wcs.all_pix2world(x_points, y_points, 1)
    x_points, y_points = coords[0], coords[1]

    x_range, y_range = max(x_points) - min(x_points), max(y_points) - min(y_points)

    plot.x_range = Range1d(max(x_points), min(x_points))
    plot.y_range = Range1d(min(y_points), max(y_points))

    """
    # true-colour colourmap
    cmap = wavelength_to_cmap(490.012)
    cmap = wavelength_to_cmap(624.127)

    plot.image(image=[image_data], x=x_points[0], y=y_points[0], dw=x_range, dh=y_range, palette=cmap, level="image", origin="bottom_right", anchor="bottom_right")
    """

    image_data = image.hdu.data

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
        wavelength = EFFECTIVE_WAVELENGTHS[image.survey][image.band]
        colour_mapper = LinearColorMapper(palette=wavelength_to_cmap(wavelength), low=vmin, high=vmax)

    if kwargs.get("relative_axes", True):
        x_arcsec_edges, y_arcsec_edges = get_relative_axes(image)
        x_range = x_arcsec_edges.max() - x_arcsec_edges.min()
        y_range = y_arcsec_edges.max() - y_arcsec_edges.min()
        focus_ra, focus_dec = 0, 0

        plot.x_range = Range1d(np.min(x_arcsec_edges), np.max(x_arcsec_edges))
        plot.y_range = Range1d(np.min(y_arcsec_edges), np.max(y_arcsec_edges))

        plot.xaxis.axis_label = "Right Ascensions / arcsec"
        plot.yaxis.axis_label = "Declination / arcsec"

        focus_ra, focus_dec = 0.0, 0.0
    else:
        focus_ra, focus_dec = image.focus.ra.value, image.focus.dec.value

    plot.image(
        image=[image_data],
        x=-x_range / 2,
        y=-y_range / 2,
        dw=x_range,
        dh=y_range,
        level="image",
        origin="bottom_left",
        anchor="bottom_left",
        color_mapper=colour_mapper,
    )

    plot.scatter(x=focus_ra, y=focus_dec, marker="cross", color="lime", size=25, line_width=4)

    return format_plot("image", plot)
