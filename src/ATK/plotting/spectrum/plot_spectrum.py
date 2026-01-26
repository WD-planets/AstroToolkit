import copy

import astropy.units as u
from bokeh.plotting import figure

from ...structures.DataSet import DataSet
from ...structures.SED import SED
from ...structures.Spectrum import Spectrum
from ..formatting import format_plot
from ..plotting_core import get_axis_label


def overlay_sed(plot: figure, spectrum: Spectrum, seds: DataSet | list[SED] | SED):
    from ..sed.plot_sed import plot as plot_sed

    if isinstance(seds, DataSet):
        seds = seds.data
    elif isinstance(seds, list):
        pass
    elif isinstance(seds, SED):
        seds = [seds]
    else:
        raise ValueError(f"Invalid type for SED overlay '{type(seds)}'.")

    for sed in seds:
        if sed._target_key != spectrum._target_key:
            continue

        sed = copy.deepcopy(sed)
        sed.flux = sed.flux.to(10**-17 * u.erg / (u.s * u.cm**2 * u.AA), equivalencies=u.spectral_density(sed.wavelength))
        sed.flux_err = sed.flux_err.to(10**-17 * u.erg / (u.s * u.cm**2 * u.AA), equivalencies=u.spectral_density(sed.wavelength))

        plot = plot_sed(sed, spectrum_plot=plot)

    return plot


def plot(spectrum: Spectrum, *args: any, **kwargs: any):
    """
    Plots an ATK Spectrum object
    """

    if not kwargs.get("sed_plot"):
        plot = figure(
            width=400,
            height=400,
            title=f"{spectrum.survey} Spectrum",
            x_axis_label=get_axis_label(spectrum, "wavelength"),
            y_axis_label=get_axis_label(spectrum, "flux"),
            tools=("pan,wheel_zoom,box_zoom,reset"),
        )
    else:
        # if plotting as an SED overlay
        plot = kwargs["sed_plot"]

    # plot spectrum
    plot.line(spectrum._get_attr_value("wavelength"), spectrum._get_attr_value("flux"), color="black", line_width=1)

    if kwargs.get("overlay"):
        plot = overlay_sed(plot, spectrum, kwargs["overlay"])

    return format_plot("spectrum", plot)
