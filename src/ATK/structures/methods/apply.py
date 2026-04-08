import numpy as np
from bokeh.models import Column, DataTable, Row
from bokeh.plotting import figure

from ..DataSet import DataSet


def unpack_layout(layout):
    plots = []

    if isinstance(layout, (figure, DataTable)):
        plots.append(layout)

    elif hasattr(layout, "children"):
        for child in layout.children:
            plots.extend(unpack_layout(child))

    return plots


def apply_methods(struct: DataSet, method: str, *args, **kwargs):
    if not struct.data:
        return struct

    data_methods = getattr(struct.data[0], "_data_methods", [])
    data_group_methods = getattr(struct.data[0], "_group_data_methods", [])
    plot_methods = getattr(struct.data[0], "_plot_methods", [])
    plot_group_methods = getattr(struct.data[0], "_group_plot_methods", [])

    all_methods = []
    for arr in [data_methods, data_group_methods, plot_methods, plot_group_methods]:
        if isinstance(arr, (list, tuple)):
            all_methods += arr
        elif isinstance(arr, dict):
            all_methods += list(arr.keys())

    if method not in all_methods:
        raise ValueError(f"Unknown plot or data method '{method}'.")

    if method in data_methods:
        data = []
        for ctnr in struct.data:
            if hasattr(ctnr, method):
                returned_ctnr = getattr(ctnr, method)(*args, **kwargs)
                data.append(returned_ctnr)
            else:
                raise ValueError(f"{type(ctnr).__name__} data does not support the method '{method}'.")
        struct.data = data

        return struct

    if method in data_group_methods:
        # collect by target key
        data = []
        keys = list(set([ctnr._target_key for ctnr in struct.data]))
        surveys = list(set([ctnr.survey for ctnr in struct.data]))

        for survey in surveys:
            survey_containers = [ctnr for ctnr in struct.data if ctnr.survey == survey]
            for key in keys:
                ctnrs = [ctnr for ctnr in survey_containers if ctnr._target_key == key]
                if not ctnrs:
                    continue
                returned_ctnrs = data_group_methods[method](ctnrs, *args, **kwargs)
                if isinstance(returned_ctnrs, list):
                    data += returned_ctnrs
                else:
                    data.append(returned_ctnrs)
        struct.data = data

        return struct

    # can't open same figure twice, so replot original data with previous plotting parameters and modify this
    figure_copy = struct.plot(**struct._stored_plot_params).figure
    if not figure:
        return struct

    plots = unpack_layout(figure_copy)
    all_figures = []

    if method in plot_methods:
        for plot in plots:
            ctnrs = [ctnr for ctnr in struct.data if ctnr._plot_id == plot.id]

            if struct._plot_method == "individual":
                # should only ever be a 1:1 mapping for ctnr:plot
                figures = [plot_methods[method](plot, ctnrs[0], **kwargs) for ctnr in ctnrs]
                all_figures.extend(figures)

            elif struct._plot_method == "combined":
                figures = plot_methods[method](plot, ctnrs, **kwargs)
                all_figures.extend(figures)

    # get rid of any None figures (shouldn't ever happen)
    all_figures = [f for f in all_figures if f is not None]

    figs_per_col = int(np.ceil(np.sqrt(len(all_figures))))

    # combines multiple plots into a grid layout of FIGS_PER_COLUMN rows and any number of columns
    all_figures = [all_figures[i : i + figs_per_col] for i in range(0, len(all_figures), figs_per_col)]
    rows = [Column(*col) for col in all_figures]

    rows = [Column(*col) for col in all_figures]
    struct.figure = Row(*rows)

    return struct
