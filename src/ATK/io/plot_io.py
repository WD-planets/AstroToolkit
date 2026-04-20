import glob
import importlib
import os
import tempfile
import time
import warnings
from pathlib import Path
from types import FunctionType

import numpy as np
from bokeh.io import output_file
from bokeh.io import save as bokeh_save
from bokeh.models import Column, Row
from bokeh.plotting import figure, show

from ..configuration.base_config import BASE_CONFIG
from ..structures.DataSet import DataSet
from ..structures.structures_core import Container
from ..structures.Target import Target
from ..utilities.mapping import build_map


def do_plotting(
    all_figures: list, plotting_func: FunctionType, plot_method: str, containers: list[Container], target: Target | None = None, **kwargs
):
    # plot .data containers individually (e.g. images)
    if plot_method == "individual":
        figures = [plotting_func(ctnr, **kwargs) for ctnr in containers]
        # flatten list if e.g. SED plotting with spectral overlay returned multiple SED plots due to having to overlay multiple spectra
        try:
            figures = [fig for fig_list in figures for fig in fig_list]
        except TypeError:
            pass

    # combine multiple .data containers into single plots (e.g. light curves)
    elif plot_method == "combined":
        figures = plotting_func(containers, **kwargs)

    # e.g. if no data was returned and plotting was attempted
    if not figures:
        return all_figures

    # add targeting info to plot titles if split by target
    if target:
        for plot in figures:
            if not hasattr(plot, "title"):
                continue
            if not plot.title:
                continue
            if target.identifier:
                plot.title.text = f"{target.identifier} {plot.title.text}"
            else:
                plot.title.text = f"{target.initial_coords.ra.value:.3f}° {target.initial_coords.dec.value:.3f} {plot.title.text}°"

    all_figures.extend(figures)

    return all_figures


def dispatch_plotting(
    all_figures: list,
    plotting_func: FunctionType,
    survey_split: bool,
    plot_method: str,
    containers: list[Container],
    target: Target | None = None,
    **kwargs,
):
    if not survey_split:
        all_figures = do_plotting(all_figures, plotting_func, plot_method, containers, target, **kwargs)

    else:
        surveys = list(set(ctnr.survey for ctnr in containers))

        for survey in surveys:
            survey_containers = [ctnr for ctnr in containers if ctnr.survey == survey]
            all_figures = do_plotting(all_figures, plotting_func, plot_method, survey_containers, target, **kwargs)

    return all_figures


def pre_plotting(structure: DataSet, kwargs: dict):
    if structure.kind == "Lightcurve":
        all_bands = sorted(list(set(lc.band for lc in structure.data)))
        kwargs["all_bands"] = all_bands

    return kwargs


def plot_data(structure: DataSet, **kwargs: any) -> figure:
    """
    Plots the data stored in a DataSet, and saves it to the .figure attribute of the data structure
    """

    if not structure.data:
        warnings.warn("Structure contains no data to for plotting.")
        return None

    ctnr_kind = structure._ctnr_kind
    module = importlib.import_module(f"ATK.plotting.{ctnr_kind}")
    plot_map = build_map(module, "plot", prefix="plot_")
    plotting_func = plot_map[ctnr_kind]

    target_keys = [t._key for t in structure.targets]

    all_figures = []

    # get any additional plotting arguments
    kwargs = pre_plotting(structure, kwargs)

    # force splitting of all containers by target
    if structure._split_by_target:
        completed_keys = []

        for target, key in zip(structure.targets, target_keys):
            if key in completed_keys:
                continue

            containers = structure._fetch_by_key(key)
            all_figures = dispatch_plotting(
                all_figures, plotting_func, structure._split_by_survey, structure._plot_method, containers, target, **kwargs
            )
            completed_keys.append(key)

    # don't split by target
    else:
        all_figures = dispatch_plotting(
            all_figures, plotting_func, structure._split_by_survey, structure._plot_method, structure.data, **kwargs
        )

    if not all_figures:
        return None

    # get rid of any None figures (shouldn't ever happen)
    all_figures = [f for f in all_figures if f is not None]

    figs_per_col = int(np.ceil(np.sqrt(len(all_figures))))

    # combines multiple plots into a grid layout of FIGS_PER_COLUMN rows and any number of columns
    all_figures = [all_figures[i : i + figs_per_col] for i in range(0, len(all_figures), figs_per_col)]
    rows = [Column(*col) for col in all_figures]

    return Row(*rows)


def check_plotted_keys(structure: DataSet, keys: list):
    if not structure._plotted_keys:
        return False

    prev_keys = sorted(structure._plotted_keys)
    if prev_keys != keys:
        return False
    else:
        return True


def open(structure: DataSet, keys: list, fname=Path | str | None, **kwargs: dict):
    """
    Opens the Bokeh plot in the .figure attribute of a DataSet in the default browser
    """

    # get previous plot parameters if they exist
    # this only really matters if a plot doesn't exist but previously did (e.g. due to using inplace=False in data methods which cannot copy a bokeh figure)
    if not kwargs and structure._stored_plot_params:
        kwargs = structure._stored_plot_params

    # determine whether to replot
    if not structure.figure:
        replot = True
    elif not check_plotted_keys(structure, keys):
        replot = True
    elif kwargs != structure._stored_plot_params:
        replot = True
    else:
        replot = False

    # plot data if it hasn't been plotted
    if replot:
        structure.plot(**kwargs)

    # if no data for plotting
    if not structure.figure:
        return

    tmp_dir = os.path.expanduser("~/.AstroToolkit/cached_figures")

    # unless a file name was provided, save to the cached_figures directory
    if not fname:
        os.makedirs(tmp_dir, exist_ok=True)

        fname_prefix = f"{structure.kind}_"

        with tempfile.NamedTemporaryFile(suffix=".html", prefix=fname_prefix, dir=tmp_dir, delete=False) as tmpfile:
            tmp_html = tmpfile.name

        output_file(tmp_html, title=structure._title)
    else:
        output_file(str(fname), title=structure._title)

    # clear cache directory of any old figures
    for f in glob.glob(os.path.join(tmp_dir, "*.html")):
        if time.time() - os.path.getmtime(f) > BASE_CONFIG._get("plot_settings", "cache_time"):
            try:
                os.remove(f)
            except PermissionError:
                pass

    show(structure.figure)


def open_basic(plot: figure, prefix: str, title: str, fname=Path | str | None):
    tmp_dir = os.path.expanduser("~/.AstroToolkit/cached_figures")

    if not fname:
        os.makedirs(tmp_dir, exist_ok=True)

        with tempfile.NamedTemporaryFile(suffix=".html", prefix=prefix, dir=tmp_dir, delete=False) as tmpfile:
            tmp_html = tmpfile.name

            output_file(tmp_html, title=title)
    else:
        output_file(fname, title=title)

    # clear cache directory of any old figures
    for f in glob.glob(os.path.join(tmp_dir, "*.html")):
        if time.time() - os.path.getmtime(f) > BASE_CONFIG._get("plot_settings", "cache_time"):
            try:
                os.remove(f)
            except PermissionError:
                pass

    show(plot)


def save(structure: DataSet, keys: list, fname=Path | str | None, **kwargs: dict):
    if not kwargs and structure._stored_plot_params:
        kwargs = structure._stored_plot_params

    # determine whether to replot
    if not structure.figure:
        replot = True
    elif not check_plotted_keys(structure, keys):
        replot = True
    elif kwargs != structure._stored_plot_params:
        replot = True
    else:
        replot = False

    # plot data if it hasn't been plotted
    if replot:
        structure.plot(**kwargs)

    # if no data for plotting
    if not structure.figure:
        return

    output_file(fname)
    bokeh_save(structure.figure, title=structure._title)
