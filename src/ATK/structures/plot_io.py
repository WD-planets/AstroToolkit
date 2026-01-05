import glob
import importlib
import os
import tempfile
import time
from pathlib import Path

from bokeh.io import output_file
from bokeh.models import Column, Row
from bokeh.plotting import figure, show

from ..configuration.base_config import BASE_CONFIG
from ..utilities.mapping import build_map
from .definitions import PlottableQueryResult

FIGS_PER_COLUMN = 3


def plot_data(kind: str, structure: PlottableQueryResult, **kwargs: any) -> figure:
    """
    Plots the data stored in a PlottableQueryResult, and saves it to the .figure attribute of the data structure
    """

    module = importlib.import_module(f"ATK.plotting.{structure.kind}")
    plot_map = build_map(module, "plot", prefix="plot_")
    kind = kind or structure.kind
    plotting_func = plot_map[kind]

    # plot .data containers individually (e.g. images)
    if structure._plot_method == "individual":
        figures = [plotting_func(ctnr, **kwargs) for ctnr in structure.data]

    # combine multiple .data containers into single plots (e.g. light curves)
    elif structure._plot_method == "combined":
        figures = plotting_func(structure.data, **kwargs)

    # e.g. if no data was returned and plotting was attempted
    if not figures:
        return None

    # get rid of any None figures (shouldn't ever happen)
    figures = [f for f in figures if f is not None]

    # combines multiple plots into a grid layout of FIGS_PER_COLUMN rows and any number of columns
    figures = [figures[i : i + FIGS_PER_COLUMN] for i in range(0, len(figures), FIGS_PER_COLUMN)]
    rows = [Column(*col) for col in figures]

    return Row(*rows)


def open(structure: PlottableQueryResult, fname=Path | str | None):
    """
    Opens the Bokeh plot in the .figure attribute of a PlottableQueryResult in the default browser
    """

    # plot data if it hasn't been plotted
    if not structure.figure:
        structure.plot()

    # unless a file name was provided, save to the cached_figures directory
    if not fname:
        tmp_dir = os.path.expanduser("~/.AstroToolkit/cached_figures")
        os.makedirs(tmp_dir, exist_ok=True)

        with tempfile.NamedTemporaryFile(
            suffix=".html", prefix=f"{structure.survey}_{structure.kind}_", dir=tmp_dir, delete=False
        ) as tmpfile:
            tmp_html = tmpfile.name

        output_file(tmp_html)
    else:
        output_file(fname)

    # clear cache directory of any old figures
    for f in glob.glob(os.path.join(tmp_dir, "*.html")):
        if time.time() - os.path.getmtime(f) > BASE_CONFIG.get("plot_settings", "cache_time"):
            os.remove(f)

    show(structure.figure)
