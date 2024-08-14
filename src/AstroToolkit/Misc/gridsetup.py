def format_grid_plots(dimensions, plots, grid_size):
    from bokeh.models import Label, Range1d
    from bokeh.models.layouts import Column as bokeh_column
    from bokeh.models.layouts import Row as bokeh_row
    from bokeh.models.widgets.tables import DataTable as bokeh_datatable
    from bokeh.plotting import figure
    from bokeh.plotting._figure import figure as bokeh_figure

    returned_plots = {}
    for plot_info in plots:
        # Bokeh figure
        if plot_info["figure"]:
            if isinstance(
                plot_info["figure"],
                (bokeh_figure, bokeh_datatable, bokeh_row, bokeh_column),
            ):
                panel = plot_info["figure"]
                panel.height, panel.width = (
                    plot_info["height"] * grid_size,
                    plot_info["width"] * grid_size,
                )
                panel.sizing_mode = "fixed"
            # ATK Structure
            else:
                if plot_info["figure"].figure:
                    panel = plot_info["figure"].figure

                    panel.height, panel.width = (
                        plot_info["height"] * grid_size,
                        plot_info["width"] * grid_size,
                    )
                    panel.sizing_mode = "fixed"
                else:
                    panel = figure(
                        height=plot_info["height"] * grid_size,
                        width=plot_info["width"] * grid_size,
                        x_axis_label="placeholder x",
                        y_axis_label="placeholder y",
                    )
                    (
                        panel.sizing_mode,
                        panel.x_range,
                        panel.y_range,
                        panel.xgrid.grid_line_color,
                        panel.ygrid.grid_line_color,
                        panel.outline_line_color,
                        panel.toolbar.logo,
                        panel.toolbar_location,
                        panel.xaxis.major_label_text_color,
                        panel.yaxis.major_label_text_color,
                        panel.xaxis.axis_label_text_color,
                        panel.yaxis.axis_label_text_color,
                    ) = (
                        "fixed",
                        Range1d(0, 10),
                        Range1d(0, 10),
                        None,
                        None,
                        None,
                        None,
                        None,
                        "white",
                        "white",
                        "white",
                        "white",
                    )

                    missing_plot_renderer = Label(
                        x=5,
                        y=5,
                        text=f"Missing {plot_info['figure'].kind} data",
                        text_align="center",
                        text_font_size="30px",
                    )
                    panel.add_layout(missing_plot_renderer)
        # None
        else:
            panel = figure(
                width=plot_info["width"] * grid_size,
                height=plot_info["height"] * grid_size,
            )
            (
                panel.outline_line_color,
                panel.toolbar.logo,
                panel.toolbar_location,
                panel.sizing_mode,
            ) = (
                None,
                None,
                None,
                "fixed",
            )

        returned_plots[plot_info["name"]] = panel
    return returned_plots
