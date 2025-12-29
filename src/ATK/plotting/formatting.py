from bokeh import events
from bokeh.models import CustomJS
from bokeh.plotting import figure

from ..configuration.base_config import BASE_CONFIG
from ..utilities.defaults import PLOT_DIMENSIONS


def format_plot(kind: str, plot: figure):
    """
    Formats the aesthetics of plots using the ATK config
    """

    # resize
    dimensions = PLOT_DIMENSIONS[kind]
    plot.width = int(BASE_CONFIG.get("plot_settings", "size")) * dimensions[0]
    plot.height = int(BASE_CONFIG.get("plot_settings", "size")) * dimensions[1]

    # output backend
    output_backend = str(BASE_CONFIG.get("plot_settings", "backend"))
    if output_backend not in ["canvas", "svg", "webgl"]:
        raise ValueError(f"Unsupported output backend {output_backend}. Accepted: canvas, svg, webgl")
    plot.output_backend = str(output_backend)

    # font sizes
    text_size = str(BASE_CONFIG.get("plot_settings", "font_size"))
    if not text_size.endswith("pt"):
        text_size += "pt"
    text_font = str(BASE_CONFIG.get("plot_settings", "font"))
    plot.axis.axis_label_text_font_size = text_size
    plot.axis.major_label_text_font_size = text_size
    if plot.title:
        plot.title.text_font_size = text_size
        plot.title.text_font = text_font
    if plot.legend:
        plot.legend.label_text_font = text_font
        plot.legend.label_text_font_size = text_size
    plot.axis.axis_label_text_font_style = "normal"
    plot.axis.axis_label_text_font = text_font
    plot.axis.major_label_text_font = text_font

    # toolbar/grid/titles
    if not BASE_CONFIG.get("plot_settings", "toolbars"):
        plot.toolbar_location = None
    if not BASE_CONFIG.get("plot_settings", "grids"):
        plot.grid.grid_line_color = None
    if not BASE_CONFIG.get("plot_settings", "titles"):
        plot.title = None

    # legend + interactivity
    if plot.legend:
        plot.legend.click_policy = "hide"

        toggle_legend_js = CustomJS(
            args=dict(leg=plot.legend[0]),
            code="""
                if (leg.visible) {
                    leg.visible = false
                    }
                else {
                    leg.visible = true
                }
        """,
        )

        plot.js_on_event(events.DoubleTap, toggle_legend_js)

    return plot
