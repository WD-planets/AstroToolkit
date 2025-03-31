class Dimensions(object):
    def __init__(self, height, width):
        self.height = height
        self.width = width


newline = "\n"


def get_plot(struct, **kwargs):
    plot_success = True
    if not struct.data:
        print("Note: Plot() missing data, suggests that no data was returned")
        struct.figure = None
        plot_success = False

    if plot_success:
        if struct.kind == "lightcurve":
            if kwargs["kind"] == "lightcurve":
                from .lightcurveplotting import plot_lightcurve

                plot = plot_lightcurve(struct, kwargs["colours"], kwargs["bands"], kwargs["timeformat"])
                dimensions = Dimensions(height=1, width=2)
            elif kwargs["kind"] == "powspec":
                # import os
                # from pathlib import Path

                from ..Timeseries.lomb_scargle import lomb_scargle

                # from ..Timeseries.pyaov import pyaov
                # from ..Timeseries.pyaov.pyaov_interface import get_analysis

                if kwargs["method"] == "ls":
                    plot, fpeak = lomb_scargle(
                        struct.data,
                        survey=struct.survey,
                        start_freq=kwargs["start_freq"],
                        stop_freq=kwargs["stop_freq"],
                        samples=kwargs["samples"],
                    ).powspec_plot
                    struct.fpeak = fpeak.value
                else:
                    raise Exception("Invalid timeseries analysis method. Accepted methods: [ls].")

                    """
                    path = Path(pyaov.__file__).parent.absolute()
                    if str(path) not in os.environ["PATH"]:
                        os.environ["PATH"] += str(path)
                    plot, fpeak = get_analysis(
                        struct=struct,
                        method=kwargs["method"],
                        gui=False,
                        start_freq=kwargs["start_freq"],
                        stop_freq=kwargs["stop_freq"],
                        samples=kwargs["samples"],
                    )
                    struct.fpeak = fpeak
                    """

                dimensions = Dimensions(height=1, width=1)

                from ..FileHandling.file_naming import generate_plotname

                generate_plotname(struct, "ATKpowspec")
            elif kwargs["kind"] == "phasefold":
                from ..Timeseries.lomb_scargle import lomb_scargle

                plot = lomb_scargle(
                    struct.data,
                    kwargs["freq"],
                    kwargs["bins"],
                    foverlay=kwargs["foverlay"],
                    repeat=kwargs["repeat"],
                    shift=kwargs["shift"],
                ).phasefold_plot
                dimensions = Dimensions(height=1, width=1)

                from ..FileHandling.file_naming import generate_plotname

                generate_plotname(struct, "ATKphasefold")

        elif struct.kind == "image":
            from .imageplotting import plot_image

            plot = plot_image(struct, kwargs["searchradius"])
            dimensions = Dimensions(height=2, width=2)
        elif struct.kind == "sed":
            from .sedplotting import plot_sed

            plot = plot_sed(struct, kwargs["spectrum_overlay"])
            dimensions = Dimensions(height=1, width=2)
        elif struct.kind == "spectrum":
            from .spectrumplotting import plot_spectrum

            plot = plot_spectrum(struct)
            dimensions = Dimensions(height=1, width=2)
        elif struct.kind == "hrd":
            from .hrdplotting import plot_hrd

            plot = plot_hrd(struct)
            dimensions = Dimensions(height=1, width=1)

        from ..Configuration.baseconfig import ConfigStruct

        if not plot:
            plot_success = False

        if plot_success:
            config = ConfigStruct()
            config.read_config()
            plot.width = int(config.unit_size) * dimensions.width
            plot.height = int(config.unit_size) * dimensions.height
            output_backend = str(config.output_backend)
            if output_backend not in ["canvas", "svg", "webgl"]:
                raise ValueError(f"Unsupported output backend {output_backend}. Accepted: canvas, svg, webgl")
            plot.output_backend = str(config.output_backend)

            text_size = str(config.font_size)
            if not text_size.endswith("pt"):
                text_size += "pt"

            text_font = str(config.font)

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

            if not config.show_toolbars:
                plot.toolbar_location = None
            if not config.show_grids:
                plot.grid.grid_line_color = None
            if not config.show_titles:
                plot.title = None

    fname = struct.plotname

    if plot_success:
        from bokeh.plotting import output_file

        output_file(f"{fname}.html")
        struct.figure = plot

    return struct
