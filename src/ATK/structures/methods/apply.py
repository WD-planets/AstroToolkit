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

    # clear figure
    struct.figure = None

    data_methods = getattr(struct.data[0], "_data_methods", [])
    data_group_methods = getattr(struct.data[0], "_group_data_methods", [])

    all_methods = []
    for arr in [data_methods, data_group_methods]:
        if isinstance(arr, (list, tuple)):
            all_methods += arr
        elif isinstance(arr, dict):
            all_methods += list(arr.keys())

    if method not in all_methods:
        raise ValueError(f"Unknown data method '{method}'.")

    if method in data_methods:
        data = []
        for ctnr in struct.data:
            if hasattr(ctnr, method):
                returned_ctnr = getattr(ctnr, method)(*args, **kwargs)
                data.append(returned_ctnr)
            else:
                raise ValueError(f"{type(ctnr).__name__} does not support the method '{method}'.")
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
