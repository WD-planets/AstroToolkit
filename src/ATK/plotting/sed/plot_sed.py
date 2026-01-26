import copy

import astropy.units as u
from bokeh.models import BasicTickFormatter, ColumnDataSource, HoverTool
from bokeh.plotting import figure

from ...structures.DataSet import DataSet
from ...structures.SED import SED
from ...structures.Spectrum import Spectrum
from ..colours import get_palette
from ..formatting import format_plot


def overlay_spectrum(sed: SED, spectra: DataSet | list[Spectrum] | Spectrum, **kwargs):
    from ..spectrum.plot_spectrum import plot as plot_spectrum

    if isinstance(spectra, DataSet):
        spectra = spectra.data
    elif isinstance(spectra, list):
        pass
    elif isinstance(spectra, SED):
        spectra = [spectra]
    else:
        raise ValueError(f"Invalid type for SED overlay '{type(spectra)}'.")

    plots = []
    for spectrum in spectra:
        if spectrum._target_key != sed._target_key:
            continue

        plot = plot_sed(sed, **kwargs)
        spectrum = copy.deepcopy(spectrum)
        spectrum.flux = spectrum.flux.to(u.mJy, equivalencies=u.spectral_density(spectrum.wavelength))
        # sed.flux_err = spectrum.flux_err.to(u.mJy, equivalencies=u.spectral_density(spectrum.wavelength))

        plots.append(plot_spectrum(spectrum, sed_plot=plot))

    return plots


def plot_sed(sed: SED, **kwargs: any):
    if kwargs.get("spectrum_plot"):
        plot = kwargs["spectrum_plot"]
    else:
        plot = figure(
            width=400,
            height=400,
            title="Spectral Energy Distribution",
            x_axis_label="Effective Wavelength / \u212b",
            y_axis_label=r"\[\text{flux / mJy}\]",
            x_axis_type="log",
            y_axis_type="log",
            tools=("pan,wheel_zoom,box_zoom,reset"),
        )

    # make ticks more readable
    plot.yaxis.formatter = BasicTickFormatter(use_scientific=False)
    plot.xaxis.major_label_overrides = {100000: r"\[10^5\]", 200000: r"\[2\times10^5\]"}
    plot.yaxis.ticker.desired_num_ticks = 5
    plot.xaxis.ticker.desired_num_ticks = 3

    # set up HoverTool
    hvr = HoverTool(
        tooltips=[
            ("survey", "@survey"),
            ("band", "@band"),
            ("wavelength", "@wavelength \u212b"),
            ("flux", f"@flux {sed.flux.unit.to_string('unicode')}"),
            ("error", f"@flux_err {sed.flux.unit.to_string('unicode')}"),
            # "/' can be impossible to differentiate if text is small
            ("separation", f"@separation {sed.separation.unit.to_string('fits')}"),
        ]
    )
    hvr.renderers = []

    data = sed.to_dataframe()

    # set up legend label column with upper limits where relevant
    data["label"] = data["survey"].astype(str)
    data.loc[data["flux_err"].isna(), "label"] += " (upper limit)"

    # get colour map (one colour for each survey)
    surveys = sorted(data["survey"].unique())
    colours = get_palette(len(surveys), shift=1)
    # get dict of mag_name: colour and map to data dataframe
    colour_map = dict(zip(surveys, colours))
    data["colour"] = data["survey"].map(colour_map)

    # split into full detections and upper limits (i.e. flux_err is nan)
    mask = data["flux_err"].isna()
    non_nan_err = data.loc[~mask].copy()
    nan_err = data.loc[mask].copy()

    # plot full detections
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

        # plot errors
        err_xs = [[x, x] for x in group["wavelength"]]
        err_ys = [[y - y_err, y + y_err] for y, y_err in zip(group["flux"], group["flux_err"])]
        plot.multi_line(err_xs, err_ys, color=colour, legend_label=label, line_width=0.5, line_cap="square")

    # plot upper limits
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

    return plot


def plot(sed: SED, **kwargs: any):
    """
    Plots an ATK SED object
    """

    if kwargs.get("overlay"):
        plots = overlay_spectrum(sed, kwargs["overlay"], **kwargs)
        return [format_plot("sed", plot) for plot in plots]

    # if plotting as an overlay for spectra
    plot = plot_sed(sed, **kwargs)

    return format_plot("sed", plot)
