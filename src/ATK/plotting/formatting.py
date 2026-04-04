from bokeh import events
from bokeh.models import CustomJS
from bokeh.plotting import figure

from ..configuration.base_config import BASE_CONFIG
from ..utilities.defaults import PLOT_DIMENSIONS


def compute_borders(text_size: str, max_tick_chars: int = 9, tick_length: int = 6, tick_standoff: int = 5, axis_standoff: int = 5):
    text_size = int(text_size[:-2])

    # Average character width approximation
    char_width = 0.6 * text_size

    # Left border estimate (y-axis footprint)
    tick_width = max_tick_chars * char_width

    left = tick_width + tick_length + tick_standoff + text_size + axis_standoff

    # Bottom border estimate (x-axis footprint)
    bottom = (3 * text_size) + tick_length + tick_standoff + text_size + axis_standoff

    return int(left), int(bottom)


def set_plot_size(plot: figure, frame_width: int, frame_height: int, text_size: str, force_square: bool):
    left, bottom = compute_borders(text_size)
    right, top = 0, int(text_size[:-2])

    plot.sizing_mode = "fixed"
    plot.width = int(frame_width) + 10
    plot.height = int(frame_height)

    if force_square:
        fw = frame_width - left - right
        fh = frame_height - top - bottom
        frame_unit = min([fw, fh])

        plot.frame_width = int(frame_unit)
        plot.frame_height = int(frame_unit)

        plot.min_border_left = left
        plot.min_border_bottom = left
        plot.min_border_top = 0
        plot.min_border_right = 10

    else:
        plot.frame_width = int(frame_width - left - right)
        plot.frame_height = int(frame_height - top - bottom)

        plot.min_border_left = left
        plot.min_border_bottom = bottom
        plot.min_border_top = top
        plot.min_border_right = right

    plot.toolbar.logo = None

    return plot


def format_plot(kind: str, plot: figure, force_square: bool = False):
    """
    Formats the aesthetics of plots using the ATK config
    """

    # resize
    dimensions = PLOT_DIMENSIONS[kind]
    width = int(BASE_CONFIG._get("plot_settings", "size")) * dimensions[0]
    height = int(BASE_CONFIG._get("plot_settings", "size")) * dimensions[1]

    # output backend
    output_backend = str(BASE_CONFIG._get("plot_settings", "backend"))
    if output_backend not in ["canvas", "svg", "webgl"]:
        raise ValueError(f"Unsupported output backend {output_backend}. Accepted: canvas, svg, webgl")
    plot.output_backend = str(output_backend)

    # font sizes
    text_size = str(BASE_CONFIG._get("plot_settings", "font_size"))
    if not text_size.endswith("pt"):
        text_size += "pt"

    plot = set_plot_size(plot, width, height, text_size, force_square)

    text_font = str(BASE_CONFIG._get("plot_settings", "font"))
    plot.axis.axis_label_text_font_size = text_size
    plot.axis.major_label_text_font_size = text_size
    if plot.title:
        plot.title.text_font_size = text_size
        plot.title.text_font = text_font
    plot.axis.axis_label_text_font_style = "normal"
    plot.axis.axis_label_text_font = text_font
    plot.axis.major_label_text_font = text_font

    # toolbar/grid/titles
    if not BASE_CONFIG._get("plot_settings", "toolbars"):
        plot.toolbar_location = None
    if not BASE_CONFIG._get("plot_settings", "grids"):
        plot.grid.grid_line_color = None
    if not BASE_CONFIG._get("plot_settings", "titles"):
        plot.title = None

    # legend + interactivity
    if len(plot.legend) > 0 and not hasattr(plot, "_legend_toggle_attached"):
        legend = plot.legend[0]

        legend.click_policy = "hide"

        toggle_legend_js = CustomJS(
            args=dict(leg=legend),
            code="""
                leg.visible = !leg.visible
            """,
        )

        # decrease font size dynamically with number of items
        step = 4
        base_size = int(text_size[:-2])
        decrement = 2

        n_items = len(legend.items)
        steps = n_items // step
        size = base_size - steps * decrement
        font_size = max(size, 8)

        legend.label_text_font_size = f"{font_size}pt"
        plot.legend.label_text_font = text_font

        plot.js_on_event(events.DoubleTap, toggle_legend_js)
        plot._legend_toggle_attached = True

    return plot
