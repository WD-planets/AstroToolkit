newline = "\n"


def exportplot(struct, fname=None):
    from bokeh.io import export_png

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()

    if not struct.figure:
        print("Note: No plot was found. One will be generated with a default configuration.")
        struct.plot()

    if fname:
        if not fname.endswith(".png"):
            fname += ".png"
    else:
        if struct.plotname.endswith(".html"):
            fname = struct.plotname[:-5] + ".png"
        else:
            fname = struct.plotname + ".html"

    if hasattr(struct, "plot"):
        if struct.figure:
            config.read_config()
            if config.enable_notifications:
                print(f"Exporting plot to PNG: {fname}{newline}")

        try:
            export_png(struct.figure, filename=fname)
        except RuntimeError:
            raise Exception(
                "Some additional setup is required for plot exporting. Platform-specific information can be found in the Bokeh documentation: https://docs.bokeh.org/en/latest/docs/user_guide/output/export.html. Note that Selenium should already be installed, as this is a dependency of AstroToolkit."
            )

    else:
        print("Note: This data structure does not support plotting.")

    return fname
