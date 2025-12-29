from bokeh.models import BasicTickFormatter, ColumnDataSource, HoverTool
from bokeh.plotting import figure

from ...structures.definitions import SED
from ..colours import get_palette
from ..formatting import format_plot


def plot(sed: SED, *args: any, **kwargs: any):
    plot = figure(
        width=400,
        height=400,
        title="Spectral Energy Distribution",
        x_axis_label="Effective Wavelength / \u212b",
        y_axis_label=r"\[\text{flux / mJy}\]",
        x_axis_type="log",
        y_axis_type="log",
        tools=("pan,wheel_zoom,reset"),
    )

    plot.yaxis.formatter = BasicTickFormatter(use_scientific=False)
    plot.xaxis.major_label_overrides = {100000: r"\[10^5\]", 200000: r"\[2\times10^5\]"}

    plot.yaxis.ticker.desired_num_ticks = 5
    plot.xaxis.ticker.desired_num_ticks = 3

    hvr = HoverTool(
        tooltips=[
            ("survey", "@survey"),
            ("band", "@band"),
            ("wavelength", "@wavelength \u212b"),
            ("flux", "@flux mJy"),
            ("error", "@flux_err mJy"),
            ("separation", "@separation arcsec"),
        ]
    )
    hvr.renderers = []

    data = sed.to_dataframe()

    data["label"] = data["survey"].astype(str)
    data.loc[data["flux_err"].isna(), "label"] += " (upper limit)"

    surveys = sorted(data["survey"].unique())
    colours = get_palette(len(surveys), shift=1)
    # get dict of mag_name: colour and map to data dataframe
    colour_map = dict(zip(surveys, colours))
    data["colour"] = data["survey"].map(colour_map)

    mask = data["flux_err"].isna()
    non_nan_err = data.loc[~mask].copy()
    nan_err = data.loc[mask].copy()

    for (survey, label, colour), group in non_nan_err.groupby(["survey", "label", "colour"]):
        scatter = plot.scatter(
            source=ColumnDataSource(group),
            x="wavelength",
            y="flux",
            size=5,
            line_color="colour",
            fill_color="colour",
            legend_label=label,
            marker="circle",
        )
        hvr.renderers.append(scatter)

        err_xs = [[x, x] for x in group["wavelength"]]
        err_ys = [[y - y_err, y + y_err] for y, y_err in zip(group["flux"], group["flux_err"])]
        plot.multi_line(err_xs, err_ys, color=colour, legend_label=label, line_width=0.5, line_cap="square")

    for (survey, label), group in nan_err.groupby(["survey", "label"]):
        scatter = plot.scatter(
            source=ColumnDataSource(group),
            x="wavelength",
            y="flux",
            size=10,
            line_color="colour",
            line_width=2,
            fill_color=None,
            legend_label=label,
            marker="+",
        )
        hvr.renderers.append(scatter)

    plot.add_tools(hvr)

    return format_plot("sed", plot)
