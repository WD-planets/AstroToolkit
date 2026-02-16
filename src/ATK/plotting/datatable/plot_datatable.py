from bokeh.models import ColumnDataSource
from bokeh.models import DataTable as bokeh_DataTable
from bokeh.models import InlineStyleSheet, TableColumn

from ...configuration.base_config import BASE_CONFIG
from ...structures.DataTable import DataTable
from ...utilities.defaults import PLOT_DIMENSIONS


def autosize_table(table, source, font_size_pt, max_height=None):
    # convert pt → pixels
    font_px = int(float(font_size_pt.replace("pt", "")) * 1.333)

    # scale row + header height
    row_height = int(font_px * 1.6)
    header_height = int(font_px * 1.8)

    table.row_height = row_height

    # number of rows
    n_rows = len(source.data.get(next(iter(source.data)), []))

    height = header_height + n_rows * row_height

    if max_height:
        height = min(height, max_height)

    table.height = height

    return table


def plot(dt: list[DataTable], **kwargs):
    source = ColumnDataSource(dt.data)
    cols = []
    for col in dt.data:
        cols.append(TableColumn(field=col, title=col))

    table = bokeh_DataTable(source=source, columns=cols, width=400, height=400)

    dimensions = PLOT_DIMENSIONS["datatable"]
    table.width = int(BASE_CONFIG.get("plot_settings", "size")) * dimensions[0]

    text_size = str(BASE_CONFIG.get("plot_settings", "font_size"))
    if not text_size.endswith("pt"):
        text_size += "pt"
    text_font = str(BASE_CONFIG.get("plot_settings", "font"))

    style_sheet = InlineStyleSheet(
        css=f".slick-header-columns {{background-color: #e0e0e0 !important;font-family: {text_font.lower()};font-size: {int(text_size[:-2])}pt; font-weight: normal}}.slick-row {{font-size: {int(text_size[:-2]) - 1}pt; font-weight: normal}}"
    )

    table.stylesheets = [style_sheet]

    height = int(BASE_CONFIG.get("plot_settings", "size")) * dimensions[1]
    table = autosize_table(table, source, text_size, height)

    dt._plot_id = table.id

    return table
