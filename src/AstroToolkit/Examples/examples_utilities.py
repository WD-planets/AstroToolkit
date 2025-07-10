from bokeh.models import PanTool

base_width = 700
base_height = 700


def format(plot, width=None, height=None, change_size=True):
    if change_size:
        plot.width = int(width * base_width)
        plot.height = int(height * base_height)
        plot.frame_width = int(width * base_width)
        plot.frame_height = int(height * base_height)
    plot.toolbar_location = None
    plot.grid.grid_line_color = None
    plot.min_border = 0
    plot.title = None
    for tool in plot.select(PanTool):
        plot.remove_tools(tool)
    return plot
