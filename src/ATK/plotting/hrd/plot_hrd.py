from importlib.resources import files

import matplotlib.cm as cm
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

from ...structures.HRD import HRD
from ..formatting import format_plot

MAG_MAP = {"Gmag": "phot_g_mean_mag", "BPmag": "phot_bp_mean_mag", "RPmag": "phot_rp_mean_mag"}


def gaussian_density_colours(x: np.ndarray, y: np.ndarray, cmap, as_hex: bool = True):
    # Stack points for KDE
    values = np.vstack([x, y])

    # KDE over all points
    kernel = stats.gaussian_kde(values)
    density = kernel(values)

    # Normalise to [0, 1]
    density /= density.max()

    # Map density to RGBA
    rgba = cmap(density)

    if not as_hex:
        return rgba

    # Convert to hex
    return ["#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255)) for r, g, b, _ in rgba]


def get_sub_cmap(name, vmin=0.0, vmax=1.0, n=256):
    base = cm.get_cmap(name)
    colors = base(np.linspace(vmin, vmax, n))

    return LinearSegmentedColormap.from_list(f"{name}_sub_{vmin}_{vmax}", colors)


def setup_background(colour_str, abs_mag_band, **kwargs):
    # -----------------
    # Backdrop plotting
    # ------------------

    colours = colour_str.split("-")

    backdrop_file = files("ATK.plotting.hrd").joinpath("backdrop_hrd_allmags.fits")

    plot = figure(
        width=400, height=400, x_axis_label=colour_str, y_axis_label=abs_mag_band, tools=("pan,wheel_zoom,box_zoom,reset"), title="GAIA HRD"
    )

    with fits.open(backdrop_file) as f:
        bg_df = Table(f[1].data).to_pandas()

    bg_df["colour"] = bg_df[MAG_MAP[colours[0]]] - bg_df[MAG_MAP[colours[1]]]
    bg_df["abs_mag"] = bg_df[MAG_MAP[abs_mag_band]] + 5 * np.log10(bg_df["parallax"] / 1000) + 5

    cmap_names = ["Greens", "Purples", "Reds", "Blues", "YlOrBr", "PuRd"]
    types = bg_df["type"].unique()
    type_to_cmap_name = dict(zip(types, cmap_names))

    bg_df["cmap"] = bg_df["type"].map(type_to_cmap_name)

    cmaps = {name: get_sub_cmap(name, 0.4, 1.0) for name in type_to_cmap_name.values()}
    bg_df["plot_colour"] = None

    for obj_type, group in bg_df.groupby("type"):
        x = group["colour"].to_numpy()
        y = group["abs_mag"].to_numpy()
        cmap = cmaps[group["cmap"].iloc[0]]

        colours = gaussian_density_colours(x, y, cmap)

        bg_df.loc[group.index, "plot_colour"] = colours
        bg_df.loc[group.index, "legend_label"] = obj_type

    bg_df["source_id"] = bg_df["source_id"].astype(str)

    if kwargs.get("background"):
        frac = np.clip(kwargs["background"], 0.0, 1.0)
        trimmed_groups = [g.sample(frac=frac, random_state=1) for _, g in bg_df.groupby("type")]
        bg_df = pd.concat(trimmed_groups, ignore_index=True)

    obj_types = {"100pc": "Stars < 100pc", "CV": "CV", "SD": "Subdwarf", "WD": "White Dwarf", "WD+dM": "WD + M-dwarf"}

    source = ColumnDataSource(bg_df)
    for obj_type, group in bg_df.groupby("type"):
        if obj_type == "ELM":
            continue

        source = ColumnDataSource(group)
        plot.scatter(
            x="colour",
            y="abs_mag",
            size=8,
            color="plot_colour",
            line_color=None,
            alpha=0.8,
            source=source,
            marker="circle",
            legend_label=obj_types[obj_type],
        )

    return plot


def overlay_source(plot: figure, hrd: HRD):
    df = pd.DataFrame({"identifier": str(hrd.identifier), "colour": hrd.colour, "abs_mag": hrd.abs_mag})
    source = ColumnDataSource(df)
    scatter = plot.scatter(
        x="colour",
        y="abs_mag",
        source=source,
        marker="square_dot",
        line_color="black",
        fill_color=None,
        size=20,
        line_width=3,
        legend_label="Gaia source(s)",
        level="overlay",
    )
    plot.y_range.flipped = True

    return plot, scatter


def plot(hrds: list[HRD], **kwargs):
    """
    Overlays multiple HRD containers on a single HRD plot
    """

    if len(set(hrd.abs_mag_band for hrd in hrds)) > 1:
        raise ValueError("Multiple HRD absolute magnitude bands detected.")
    if len(set(hrd.colour_bands for hrd in hrds)) > 1:
        raise ValueError("Multiple HRD colour bands detected.")

    abs_mag_band = hrds[0].abs_mag_band
    colours = hrds[0].colour_bands

    figs = []
    if kwargs.get("split", True):
        for _ in hrds:
            figs.append(setup_background(colours, abs_mag_band, **kwargs))
    else:
        figs.append(setup_background(colours, abs_mag_band, **kwargs))

    # --------------
    # Source Overlay
    # --------------

    if not kwargs.get("split", True):
        hvr = HoverTool(tooltips=[("id", "@identifier"), ("colour", "@colour"), ("abs_mag", "@abs_mag")])
        hvr.renderers = []
        for hrd in hrds:
            figs[0], scatter = overlay_source(figs[0], hrd)
            hvr.renderers.append(scatter)
        figs[0].add_tools(hvr)

        for hrd in hrds:
            hrd._plot_id = figs[0].id
    else:
        for fig, hrd in zip(figs, hrds):
            hvr = HoverTool(tooltips=[("id", "@identifier"), ("colour", "@colour"), ("abs_mag", "@abs_mag")])
            hvr.renderers = []
            fig, scatter = overlay_source(fig, hrd)
            hvr.renderers.append(scatter)
            fig.add_tools(hvr)
            hrd._plot_id = fig.id

    return [format_plot("hrd", plot) for plot in figs]
