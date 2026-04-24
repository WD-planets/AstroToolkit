from bokeh.models import PanTool
from bokeh.plotting import figure

GRID_SIZE = 200
TEXT_SIZE = "8pt"


def compute_borders(max_tick_chars: int = 9, tick_length: int = 6, tick_standoff: int = 5, axis_standoff: int = 5):
    text_size = int(TEXT_SIZE[:-2])

    # Average character width approximation
    char_width = 0.6 * text_size

    # Left border estimate (y-axis footprint)
    tick_width = max_tick_chars * char_width

    left = tick_width + tick_length + tick_standoff + text_size + axis_standoff

    # Bottom border estimate (x-axis footprint)
    bottom = (3 * text_size) + tick_length + tick_standoff + text_size + axis_standoff

    return int(left), int(bottom)


def set_panel_size(panel: figure, force_square: bool, height: int, width: int, shift_outline: bool = False):
    left, bottom = compute_borders()
    right, top = 0, int(TEXT_SIZE[:-2])

    frame_width = GRID_SIZE * width
    frame_height = GRID_SIZE * height

    panel.sizing_mode = "fixed"
    panel.width = int(frame_width) + 10
    panel.height = int(frame_height)

    if force_square:
        fw = frame_width - left - right
        fh = frame_height - top - bottom
        frame_unit = min([fw, fh])

        panel.frame_width = int(frame_unit)
        panel.frame_height = int(frame_unit)

        panel.min_border_left = left
        panel.min_border_bottom = left
        panel.min_border_top = 0
        panel.min_border_right = 10

    elif shift_outline:
        panel.frame_width = int(frame_width - left - right)
        panel.frame_height = int(frame_height - top - bottom)

        panel.min_border_left = left
        panel.min_border_bottom = left
        panel.min_border_top = 0
        panel.min_border_right = 10

    else:
        panel.frame_width = int(frame_width - left - right)
        panel.frame_height = int(frame_height - top - bottom)

        panel.min_border_left = left
        panel.min_border_bottom = bottom
        panel.min_border_top = top
        panel.min_border_right = right

    panel.toolbar.logo = None

    return panel


def format_plot(plot, width: float, height: float, force_square: bool = False):
    """
    Recursively formats a Bokeh figure/layout
    """

    def _format_single(fig):
        if isinstance(fig, figure):
            fig.sizing_mode = "fixed"

            fig.toolbar_location = None
            fig.grid.grid_line_color = None
            fig.title = None

            fig = set_panel_size(fig, force_square, height, width)

            for tool in fig.select(PanTool):
                fig.remove_tools(tool)

    def _recurse_layout(obj):
        if isinstance(obj, (figure)):
            _format_single(obj)
            return obj

        elif hasattr(obj, "children"):
            for child in obj.children:
                _recurse_layout(child)
            return obj

        else:
            return obj

    return _recurse_layout(plot)
