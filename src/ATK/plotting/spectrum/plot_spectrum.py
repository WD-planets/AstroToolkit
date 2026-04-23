import copy

import astropy.units as u
from bokeh.plotting import figure

from ...structures.DataSet import DataSet
from ...structures.methods.spectrum.fitting import do_fitting
from ...structures.methods.spectrum.radial_velocities import get_rvs
from ...structures.SED import SED
from ...structures.Spectrum import Spectrum
from ..formatting import format_plot
from ..plotting_core import get_axis_label
from .spectrum_overlay import plot_overlay


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

    if spectrum.wavelength is not None:
        x = "wavelength"
    elif spectrum.velocity is not None:
        x = "velocity"
    else:
        raise ValueError("Spectrum requires one of 'wavelength', 'velocity'.")

    if not kwargs.get("sed_plot"):
        plot = figure(
            width=400,
            height=400,
            title=f"{spectrum.survey.upper()} Spectrum",
            x_axis_label=get_axis_label(spectrum, x),
            y_axis_label=get_axis_label(spectrum, "flux"),
            tools=("pan,wheel_zoom,box_zoom,reset"),
        )
    else:
        # if plotting as an SED overlay
        plot = kwargs["sed_plot"]

    # plot spectrum
    plot.line(spectrum._get_attr_value(x), spectrum._get_attr_value("flux"), color="black", line_width=1)

    if kwargs.get("overlay"):
        plot = overlay_sed(plot, spectrum, kwargs["overlay"])

    plot = plot_overlay(plot, spectrum)

    spectrum._plot_id = plot.id

    # Fitting
    # -------

    if kwargs.get("fit") or kwargs.get("rv_fit"):
        if spectrum.wavelength is None and spectrum.velocity is not None:
            raise Exception("Fitting is not supported for velocity spectra.")

    if kwargs.get("fit"):
        plot = do_fitting(plot, spectrum, **kwargs)

    if kwargs.get("rv_fit"):
        plot = get_rvs(plot, spectrum, **kwargs)

    return format_plot("spectrum", plot)
