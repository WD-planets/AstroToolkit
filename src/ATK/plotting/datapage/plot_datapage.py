from bokeh.models import GridBox, InlineStyleSheet, Label, Range1d
from bokeh.plotting import figure

from ...configuration.base_config import BASE_CONFIG
from ...plotting.datatable.plot_datatable import autosize_table
from ...structures.DataPages import DataPages
from ...structures.DataSet import DataSet
from ...structures.methods.apply import unpack_layout

TEXT_SIZE = str(BASE_CONFIG._get("datapage_settings", "font_size"))
if not TEXT_SIZE.endswith("pt"):
    TEXT_SIZE = f"{TEXT_SIZE}pt"
TEXT_FONT = str(BASE_CONFIG._get("datapage_settings", "font"))


def format_datatable(table, height, width):
    grid_size = BASE_CONFIG._get("datapage_settings", "grid_size")

    table.min_width = 0
    table = autosize_table(table, table.source, TEXT_SIZE, grid_size * height)

    style_sheet = InlineStyleSheet(
        css=f".slick-header-columns {{background-color: #e0e0e0 !important;font-family: {TEXT_FONT.lower()};font-size: {int(TEXT_SIZE[:-2])}pt; font-weight: normal}}.slick-row {{font-size: {int(TEXT_SIZE[:-2]) - 1}pt; font-weight: normal}}"
    )

    table.stylesheets = [style_sheet]

    return table


def compute_borders(max_tick_chars: int = 9, tick_length: int = 6, tick_standoff: int = 5, axis_standoff: int = 5):
    text_size = int(TEXT_SIZE[:-2])

    # average character width approximation
    char_width = 0.6 * text_size

    # left border estimate
    tick_width = max_tick_chars * char_width

    left = tick_width + tick_length + tick_standoff + text_size + axis_standoff

    # bottom border estimate
    bottom = (3 * text_size) + tick_length + tick_standoff + text_size + axis_standoff

    return int(left), int(bottom)


def set_font_sizes(panel):
    panel.axis.axis_label_text_font_size = TEXT_SIZE
    panel.axis.major_label_text_font_size = TEXT_SIZE
    if panel.title:
        panel.title.text_font_size = TEXT_SIZE
        panel.title.text_font = TEXT_FONT
    if panel.legend:
        panel.legend.label_text_font = TEXT_FONT
        panel.legend.label_text_font_size = TEXT_SIZE
    panel.axis.axis_label_text_font_style = "normal"
    panel.axis.axis_label_text_font = TEXT_FONT
    panel.axis.major_label_text_font = TEXT_FONT

    return panel


def set_panel_size(panel: figure, force_square: bool, height: int, width: int, shift_outline: bool = False):
    grid_size = BASE_CONFIG._get("datapage_settings", "grid_size")

    left, bottom = compute_borders()
    right, top = 0, int(TEXT_SIZE[:-2])

    frame_width = grid_size * width
    frame_height = grid_size * height

    panel.sizing_mode = "fixed"
    panel.width = frame_width
    panel.height = frame_height

    if force_square:
        fw = frame_width - left - right
        fh = frame_height - top - bottom
        frame_unit = min([fw, fh])

        panel.frame_width = frame_unit
        panel.frame_height = frame_unit

        panel.min_border_left = left
        panel.min_border_bottom = left
        panel.min_border_top = 0
        panel.min_border_right = 0

    elif shift_outline:
        panel.frame_width = frame_width - left - right
        panel.frame_height = frame_height - top - bottom

        panel.min_border_left = left
        panel.min_border_bottom = left
        panel.min_border_top = 0
        panel.min_border_right = 0

    else:
        panel.frame_width = frame_width - left - right
        panel.frame_height = frame_height - top - bottom

        panel.min_border_left = left
        panel.min_border_bottom = bottom
        panel.min_border_top = top
        panel.min_border_right = right

    panel.toolbar.logo = None

    return panel


def hide_panel_visuals(panel, outline=False):
    if not outline:
        panel.outline_line_color = None
    panel.toolbar_location = None
    panel.toolbar.logo = None
    panel.grid.grid_line_color = None
    panel.xaxis.visible = False
    panel.yaxis.visible = False
    panel.xaxis.major_label_text_color = "white"
    panel.yaxis.major_label_text_color = "white"
    panel.xaxis.axis_label_text_color = "white"
    panel.yaxis.axis_label_text_color = "white"

    return panel


def blank_panel(height, width):
    panel = figure()
    panel = hide_panel_visuals(panel)

    # suppress MISSING_RENDERERS warning
    panel.scatter(x=[], y=[], visible=False)

    return panel


def get_missing_panel(dataset, height, width):
    panel = figure(x_axis_label="placeholder x", y_axis_label="placeholder y")
    panel.x_range = Range1d(0, 10)
    panel.y_range = Range1d(0, 10)
    missing_plot_renderer = Label(x=5, y=5, text=f"Missing {dataset.kind} data", text_align="center", text_font_size="30px")
    panel.add_layout(missing_plot_renderer)
    panel = hide_panel_visuals(panel, outline=True)

    return panel


def parse_layout(grid):
    """
    Parses a 2D grid of objects (or None) into rectangular regions
    """

    nrows = len(grid)
    ncols = len(grid[0])

    visited = [[False] * ncols for _ in range(nrows)]
    regions = []

    for r in range(nrows):
        for c in range(ncols):
            obj = grid[r][c]

            if visited[r][c]:
                continue

            # Only start region if this is the top-left boundary (unless None for blank panel)
            if obj is not None and ((r > 0 and grid[r - 1][c] is obj) or (c > 0 and grid[r][c - 1] is obj)):
                continue

            # Expand horizontally
            colspan = 0
            while c + colspan < ncols and grid[r][c + colspan] is obj:
                colspan += 1

            # Expand vertically
            rowspan = 0
            while r + rowspan < nrows and all(grid[r + rowspan][c + k] is obj for k in range(colspan)):
                rowspan += 1

            # Mark as visited
            for rr in range(r, r + rowspan):
                for cc in range(c, c + colspan):
                    visited[rr][cc] = True

            regions.append(
                {
                    "object_id": id(obj) if obj is not None else None,
                    "object": obj,
                    "row": r,
                    "col": c,
                    "rowspan": rowspan,
                    "colspan": colspan,
                }
            )

    return regions


def validate_layout(layout: list[list[DataSet]]):
    regions = parse_layout(layout)
    seen = {}

    for region in regions:
        obj_id = region["object_id"]
        if obj_id in seen:
            raise ValueError("Dataset appears in multiple regions.")
        seen[obj_id] = region


def prepare_datasets(key: str, datasets: list[DataSet]):
    plot_dict = {}
    for dataset in datasets:
        ctnrs = [ctnr for ctnr in dataset.data if ctnr._target_key == key]
        if not dataset.figure:
            if dataset.kind == "hrd":
                dataset.plot(combine=False)
            else:
                dataset.plot()
        plot_ids = [ctnr._plot_id for ctnr in ctnrs]

        plots = unpack_layout(dataset.figure)
        plots = [plot for plot in plots if plot.id in plot_ids]

        # need to come up with a proper solution for this
        if len(plots) > 1:
            plots = plots[0:1]

        plot_dict[id(dataset)] = plots

    return plot_dict


def get_datapage(layout: list[list]):
    datasets = list({id(ds): ds for row in layout for ds in row if ds is not None}.values())
    target_keys = list({target._key for ds in datasets for target in ds.targets})

    validate_layout(layout)
    regions = parse_layout(layout)

    dps = []
    plot_map = {}
    for index, key in enumerate(target_keys):
        plots = prepare_datasets(key, datasets)

        grid_children = []

        for region in regions:
            obj_id = region["object_id"]

            force_square, shift_outline = False, False
            if obj_id is None:
                plot = blank_panel(region["rowspan"], region["colspan"])
                kind = None
            else:
                kind = region["object"].kind

                if not plots[obj_id]:
                    plot = get_missing_panel(region["object"], region["rowspan"], region["colspan"])
                    shift_outline = True
                else:
                    plot = plots[obj_id][0]
                    force_square = True if kind in ["image"] else False

            if kind not in ["datatable"]:
                plot = set_panel_size(plot, force_square, region["rowspan"], region["colspan"], shift_outline)
                plot = set_font_sizes(plot)
            else:
                plot = format_datatable(plot, region["rowspan"], region["colspan"])

            grid_children.append((plot, region["row"], region["col"], region["rowspan"], region["colspan"]))

        final_layout = GridBox(children=grid_children)

        plot_map[key] = final_layout.id
        dps.append(final_layout)

    targets = []
    for ds in datasets:
        existing_keys = [target._key for target in targets]
        for target in ds.targets:
            if target._key not in existing_keys:
                targets.append(target)

    datapages = DataPages(targets=targets, figures=dps, _plot_map=plot_map)

    return datapages
