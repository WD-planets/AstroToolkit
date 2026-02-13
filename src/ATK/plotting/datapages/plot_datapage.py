from bokeh.io import show
from bokeh.layouts import Column, Row

from ...structures.DataSet import DataSet
from ...structures.methods.apply import unpack_layout

BASE_WIDTH = 300
BASE_HEIGHT = 300


def parse_layout(grid):
    nrows = len(grid)
    ncols = len(grid[0])
    visited = [[False] * ncols for _ in range(nrows)]
    regions = []

    for r in range(nrows):
        for c in range(ncols):
            if visited[r][c]:
                continue

            label = grid[r][c]

            colspan = 0
            while c + colspan < ncols and grid[r][c + colspan] == label:
                colspan += 1

            rowspan = 0
            while r + rowspan < nrows and all(grid[r + rowspan][c + k] == label for k in range(colspan)):
                rowspan += 1

            for rr in range(r, r + rowspan):
                for cc in range(c, c + colspan):
                    visited[rr][cc] = True

            regions.append({"name": label, "row": r, "col": c, "rowspan": rowspan, "colspan": colspan})

    return regions


def prepare_datasets(key: str, datasets: list[DataSet]):
    plot_dict = {}
    for dataset in datasets:
        ctnrs = [ctnr for ctnr in dataset.data if ctnr._target_key == key]
        if not dataset.figure:
            dataset.plot()
        plot_ids = [ctnr._plot_id for ctnr in ctnrs]

        plots = unpack_layout(dataset.figure)
        plots = [plot for plot in plots if plot.id in plot_ids]

        # need to come up with a proper solution for this
        if len(plots) > 1:
            plots = plots[0:1]

        plot_dict[dataset.kind] = plots

    return plot_dict


def get_datapage(layout: list[list], datasets: list[DataSet]):
    target_keys = []
    for dataset in datasets:
        target_keys.extend([target._key for target in dataset.targets])
    target_keys = list(set(target_keys))

    regions = parse_layout(layout)
    for key in target_keys:
        plots = prepare_datasets(key, datasets)

        rows_dict = {}
        for region in regions:
            r = region["row"]
            rows_dict.setdefault(r, []).append(region)

        rows_list = []

        for r in sorted(rows_dict):
            row_regions = rows_dict[r]

            row_regions.sort(key=lambda x: x["col"])

            row_children = []
            for region in row_regions:
                plot = plots[region["name"]][0]
                plot.width = BASE_WIDTH * region["colspan"]
                plot.height = BASE_HEIGHT * region["rowspan"]
                row_children.append(plot)

            rows_list.append(Row(*row_children))

        final_layout = Column(*rows_list)

        show(final_layout)
