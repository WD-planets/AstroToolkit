from bokeh.io import output_file, show
from bokeh.models import Label, Range1d
from bokeh.models.layouts import Column as bokeh_column
from bokeh.models.layouts import Row as bokeh_row
from bokeh.models.widgets.tables import DataTable as bokeh_datatable
from bokeh.plotting import figure
from bokeh.plotting._figure import figure as bokeh_figure

from ..Configuration.baseconfig import ConfigStruct
from ..Data.dataquery import DataStruct
from ..Data.hrdquery import HrdStruct
from ..Data.imagequery import ImageStruct
from ..Data.lightcurvequery import LightcurveStruct
from ..Data.sedquery import SedStruct
from ..Data.spectrumquery import SpectrumStruct

data_structures = (DataStruct, HrdStruct, ImageStruct, LightcurveStruct, SedStruct, SpectrumStruct)


class Datapage(object):
    """Datapage()
    This structure is returned from the :func:`datapage() <AstroToolkit.Datapages.datapage>` tool.

    .. rubric:: Attributes
        :heading-level: 1

    figure : *bokeh layout*
        the stored datapage


    |

    """

    def __init__(self, figure):
        self.figure = figure

    def __str__(self):
        return "<ATK Datapage>"

    def showplot(self, fname=None):
        """showplot(fname)

        Opens the datapage stored in the ‘figure’ attribute in the default web browser, and saves it to local files.

        :param fname: file name to save the datapage to
        :type fname: str

        |

        """
        if not fname:
            raise ValueError("Datapage showplot() requires fname.")

        if not fname.endswith(".html"):
            fname += ".html"
        output_file(fname)
        show(self.figure)

    def saveplot(self, fname=None):
        """saveplot(fname)

        Saves the datapage stored in the 'figure' attribute to local files without opening it in the web browser.

        :param fname: file name to save the figure to
        :type fname: str

        |

        """

        if not fname:
            raise ValueError("Datapage saveplot() requires fname.")

        if not fname.endswith(".html"):
            fname += ".html"
        output_file(fname)
        show(self.figure)


config = ConfigStruct()
config.read_config()
text_font = config.font
text_size = config.datapage_font_size

if not text_size.endswith("pt"):
    text_size += "pt"

excess_scaling = (int(text_size[:-2]) - 10) / 10
horizontal_padding = 80 + int(80 * excess_scaling / 1.6)
vertical_padding = 80

size_scaling = 1
frame_scaling = 1
scaling_mode = "fixed"

if not text_size.endswith("pt"):
    text_size += "pt"

text_style = "normal"


def hide_panel_visuals(panel):
    panel.outline_colour = None
    panel.toolbar_location = None
    panel.toolbar.logo = None
    panel.grid.grid_line_color = None
    panel.xaxis.major_label_text_color = "white"
    panel.yaxis.major_label_text_color = "white"
    panel.xaxis.axis_label_text_color = "white"
    panel.yaxis.axis_label_text_color = "white"
    return panel


def set_panel_size(panel, height, width):
    panel.frame_width = width
    panel.frame_height = height
    panel.min_border_left = horizontal_padding
    panel.min_border_bottom = vertical_padding
    panel.min_border_top = 0
    panel.min_border_right = 0
    panel.toolbar.logo = None
    return panel


def get_blank_panel(height, width):
    panel = figure()
    panel.frame_height = height
    panel.frame_width = width
    panel = hide_panel_visuals(panel)
    return panel


def get_missing_panel(entry, height, width):
    panel = figure(x_axis_label="placeholder x", y_axis_label="placeholder y")
    panel = set_panel_size(panel, height, width)
    panel.x_range = Range1d(0, 10)
    panel.y_range = Range1d(0, 10)
    missing_plot_renderer = Label(
        x=5, y=5, text=f"Missing {entry['figure'].kind} data", text_align="center", text_font_size="30px"
    )
    panel.add_layout(missing_plot_renderer)
    return panel


def format_plot(entry, height, width):
    panel = entry["figure"]
    panel = set_panel_size(panel, height, width)
    panel.axis.axis_label_text_font_size = text_size
    panel.axis.major_label_text_font_size = text_size
    panel.axis.axis_label_text_font_style = text_style
    panel.axis.axis_label_text_font = text_font
    panel.axis.major_label_text_font = text_font
    if panel.title:
        panel.title.text_font_size = text_size
        panel.title.text_font = text_font
    if panel.legend:
        panel.legend.label_text_font = text_font
        panel.legend.label_text_font_size = text_size
    return panel


def setup_grid(figures, datapage_layout):
    from bokeh.layouts import column, layout, row

    finished_rows = []
    for datapage_row in datapage_layout:
        row_panels = []
        for index, entry in enumerate(datapage_row):
            if isinstance(entry, list):
                column_panels = []
                for label in entry:
                    column_panels.append(figures[label])
                row_panels.append(column(column_panels))
            else:
                row_panels.append(figures[entry])
        finished_rows.append(row(row_panels))
    datapage = layout(finished_rows)

    return datapage


def format_grid_plots(dimensions, plots, grid_size, layout):
    returned_plots = {}

    for entry in plots:
        panel_height = size_scaling * entry["height"] * grid_size
        panel_width = size_scaling * entry["width"] * grid_size
        panel_name = entry["name"]

        if isinstance(entry["figure"], data_structures):
            kind = "structure"
        else:
            kind = "normal"

        if kind == "normal":
            # figure = None
            if not entry["figure"]:
                panel = get_blank_panel(panel_height, panel_width)
            # Columns/Rows/Datatables
            elif isinstance(entry["figure"], (bokeh_column, bokeh_row, bokeh_datatable)):
                panel = entry["figure"]

                if isinstance(entry["figure"], bokeh_datatable):
                    from ..DatapageElements.metadata_table import \
                        DataTable_Style

                    panel.stylesheets = [DataTable_Style().datapage_style_sheet]
            # Bokeh figures
            elif isinstance(entry["figure"], bokeh_figure):
                panel = format_plot(entry, panel_height, panel_width)
        elif kind == "structure":
            # Missing data panels
            if not entry["figure"].figure:
                panel = get_missing_panel(figure, panel_height, panel_width)
            # Data panels
            else:
                entry["figure"] = entry["figure"].figure
                panel = format_plot(entry, panel_height, panel_width)

        # Global settings
        panel.sizing_mode = scaling_mode
        panel.height = panel_height
        panel.width = panel_width

        returned_plots[panel_name] = panel

    datapage = setup_grid(returned_plots, layout)

    return Datapage(datapage)
