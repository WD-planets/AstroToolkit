import glob
import importlib
import os
import tempfile
import time
import warnings
from pathlib import Path
from types import FunctionType

from bokeh.io import output_file
from bokeh.models import Column, Row
from bokeh.plotting import figure, show

from ..configuration.base_config import BASE_CONFIG
from ..structures.DataSet import DataSet
from ..structures.structures_core import BaseContainer
from ..structures.Target import Target
from ..utilities.mapping import build_map

FIGS_PER_COLUMN = 3


def dispatch_plotting(
    all_figures: list,
    plotting_func: FunctionType,
    structure: DataSet,
    containers: list[BaseContainer],
    target: Target | None = None,
    **kwargs,
):
    # plot .data containers individually (e.g. images)
    if structure._plot_method == "individual":
        figures = [plotting_func(ctnr, **kwargs) for ctnr in containers]
        # flatten list if e.g. SED plotting with spectral overlat returned multiple SED plots due to having to overlay multiple spectra
        try:
            figures = [fig for fig_list in figures for fig in fig_list]
        except TypeError:
            pass

    # combine multiple .data containers into single plots (e.g. light curves)
    elif structure._plot_method == "combined":
        figures = plotting_func(containers, **kwargs)

    # e.g. if no data was returned and plotting was attempted
    if not figures:
        return all_figures

    # add targeting info to plot titles if split by target
    if target:
        for plot in figures:
            if target.identifier:
                plot.title.text = f"{target.identifier} {plot.title.text}"
            else:
                plot.title.text = f"{target.initial_coords.ra.value:.3f}° {target.initial_coords.dec.value:.3f}°"

    all_figures += figures

    return all_figures


def plot_data(kind: str, structure: DataSet, **kwargs: any) -> figure:
    """
    Plots the data stored in a DataSet, and saves it to the .figure attribute of the data structure
    """

    module = importlib.import_module(f"ATK.plotting.{structure.kind}")
    plot_map = build_map(module, "plot", prefix="plot_")
    kind = kind or structure.kind
    plotting_func = plot_map[kind]

    if not structure.data:
        warnings.warn("Structure contains no data to for plotting.")
        return None

    target_keys = [t._key for t in structure.targets]

    all_figures = []

    # force splitting of all plots by targets
    if structure._split_by_target:
        for target, key in zip(structure.targets, target_keys):
            containers = structure._fetch_by_key(key)
            all_figures = dispatch_plotting(all_figures, plotting_func, structure, containers, target, **kwargs)

    # don't split by target
    else:
        all_figures = dispatch_plotting(all_figures, plotting_func, structure, structure.data, **kwargs)

    if not all_figures:
        return None

    # get rid of any None figures (shouldn't ever happen)
    all_figures = [f for f in all_figures if f is not None]

    # combines multiple plots into a grid layout of FIGS_PER_COLUMN rows and any number of columns
    all_figures = [all_figures[i : i + FIGS_PER_COLUMN] for i in range(0, len(all_figures), FIGS_PER_COLUMN)]
    rows = [Column(*col) for col in all_figures]

    return Row(*rows)


def open(structure: DataSet, fname=Path | str | None, **kwargs: dict):
    """
    Opens the Bokeh plot in the .figure attribute of a DataSet in the default browser
    """

    # get previous plot parameters if they exist
    # this only really matters if a plot doesn't exist but previously did (e.g. due to using inplace=False in data methods which cannot copy a bokeh figure)
    if not kwargs and structure._stored_plot_params:
        kwargs = structure._stored_plot_params

    # plot data if it hasn't been plotted
    if not structure.figure:
        structure.plot(**kwargs)

    # if no data for plotting
    if not structure.figure:
        return

    tmp_dir = os.path.expanduser("~/.AstroToolkit/cached_figures")

    # unless a file name was provided, save to the cached_figures directory
    if not fname:
        os.makedirs(tmp_dir, exist_ok=True)

        fname_prefix = f"{structure.survey}_{structure.kind}_" if structure.survey else f"{structure.kind}_"

        with tempfile.NamedTemporaryFile(suffix=".html", prefix=fname_prefix, dir=tmp_dir, delete=False) as tmpfile:
            tmp_html = tmpfile.name

        output_file(tmp_html, title=structure._title)
    else:
        output_file(fname, title=structure._title)

    # clear cache directory of any old figures
    for f in glob.glob(os.path.join(tmp_dir, "*.html")):
        if time.time() - os.path.getmtime(f) > BASE_CONFIG.get("plot_settings", "cache_time"):
            os.remove(f)

    show(structure.figure)
