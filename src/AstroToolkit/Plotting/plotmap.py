from bokeh.io import output_file


class Dimensions(object):
    def __init__(self, height, width):
        self.height = height
        self.width = width


newline = "\n"


def saveplot(self, fname=None):
    from bokeh.plotting import save as bokeh_save

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()

    if self.figure:
        if fname:
            if not fname.endswith(".html"):
                fname += ".html"
            output_file(fname)
        else:
            fname = self.plotname
            output_file(fname)

        if hasattr(self, "plot"):
            if self.figure:
                config.read_config()
                if config.enable_notifications == "True":
                    print(f"Saving plot to local storage: {fname}{newline}")
                bokeh_save(self.figure)
            else:
                print("Note: No plot to save.")
        else:
            print("Note: No plot to save. Create one using .plot()")
    else:
        print("Note: No plot to save. Create one using .plot()")

    return fname


def showplot(self, fname=None):
    from bokeh.plotting import show as bokeh_show

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()

    if self.figure:
        if fname:
            if not fname.endswith(".html"):
                fname += ".html"
            output_file(fname)
        else:
            fname = self.plotname
            output_file(fname)

        if hasattr(self, "plot"):
            if self.figure:
                config.read_config()
                if config.enable_notifications == "True":
                    print(f"Saving plot to local storage: {fname}{newline}")

                bokeh_show(self.figure)
            else:
                print(
                    "Note: No plot to show. Suggests necesesary data was not retrieved."
                )
        else:
            print("Note: No plot to show. Create one using .plot()")
    else:
        print("Note: No plot to show. Create one using .plot()")

    return fname


def map_to_plot(
    struct,
    kind=None,
    colours=None,
    bands=None,
    spectrum_overlay=False,
    survey=None,
    freq=None,
    bins=None,
    timeformat=None,
    method=None,
    foverlay=True,
    repeat=1,
    shift=0,
):
    plot_success = True
    if not struct.data:
        print("Note: Plot() missing data, suggests that no data was returned")
        struct.plot = None
        plot_success = False

    if plot_success:
        if struct.kind == "lightcurve":
            if kind == "lightcurve":
                from .lightcurveplotting import plot_lightcurve

                plot = plot_lightcurve(struct, colours, bands, timeformat)
                dimensions = Dimensions(height=1, width=2)
            elif kind == "powspec":
                import os
                from pathlib import Path

                from ..Timeseries.lomb_scargle import lomb_scargle
                from ..Timeseries.pyaov import pyaov
                from ..Timeseries.pyaov.pyaov_interface import get_analysis

                if method == "ls":
                    plot, fpeak = lomb_scargle(
                        struct.data, survey=struct.survey
                    ).powspec_plot
                    struct.fpeak = fpeak.value
                else:
                    path = Path(pyaov.__file__).parent.absolute()
                    if str(path) not in os.environ["PATH"]:
                        os.environ["PATH"] += str(path)
                    plot, fpeak = get_analysis(struct=struct, method=method, gui=False)
                    struct.fpeak = fpeak

                dimensions = Dimensions(height=1, width=1)

                from ..FileHandling.file_naming import generate_plotname

                generate_plotname(struct, "ATKpowspec")
            elif kind in ["phase", "phasefold", "fold"]:
                from ..Timeseries.lomb_scargle import lomb_scargle

                try:
                    plot = lomb_scargle(
                        struct.data,
                        freq,
                        bins,
                        foverlay=foverlay,
                        repeat=repeat,
                        shift=shift,
                    ).phasefold_plot
                except:
                    plot = None

                dimensions = Dimensions(height=1, width=1)

                from ..FileHandling.file_naming import generate_plotname

                generate_plotname(struct, "ATKphasefold")
            else:
                raise Exception("Invalid plotting kind for ATK Lightcurve data.")
        elif struct.kind == "image":
            from .imageplotting import plot_image

            plot = plot_image(struct)
            dimensions = Dimensions(height=2, width=2)
        elif struct.kind == "sed":
            from .sedplotting import plot_sed

            plot = plot_sed(struct, spectrum_overlay, survey)
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
                raise ValueError(
                    f"Unsupported output backend {output_backend}. Accepted: canvas, svg, webgl"
                )
            plot.output_backend = str(config.output_backend)

            text_size = str(config.font_size)
            if not text_size.endswith("pt"):
                text_size += "pt"

            text_font = str(config.font)

            plot.axis.axis_label_text_font_size = text_size
            plot.axis.major_label_text_font_size = text_size
            plot.legend.label_text_font_size = text_size
            plot.title.text_font_size = text_size
            plot.legend.label_text_font = text_font
            plot.axis.axis_label_text_font_style = "normal"
            plot.axis.axis_label_text_font = text_font
            plot.axis.major_label_text_font = text_font
            plot.title.text_font = text_font

    fname = struct.plotname

    if plot_success:
        from bokeh.plotting import output_file

        output_file(f"{fname}.html")
        struct.figure = plot

    import types

    struct.showplot = types.MethodType(showplot, struct)

    struct.saveplot = types.MethodType(saveplot, struct)

    return struct
