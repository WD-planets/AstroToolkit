from bokeh.io import output_file

newline = "\n"


def showplot(struct, fname=None):
    from bokeh.plotting import show as bokeh_show

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()

    if not struct.figure:
        print("Note: No plot was found. One will be generated with a default configuration.")
        struct.plot()

    if fname:
        if not fname.endswith(".html"):
            fname += ".html"
        output_file(fname)
    else:
        fname = struct.plotname
        if not fname.endswith(".html"):
            fname += ".html"
        output_file(fname)

    if hasattr(struct, "plot"):
        if struct.figure:
            config.read_config()
            if config.enable_notifications:
                print(f"Saving plot to local storage: {fname}{newline}")

            bokeh_show(struct.figure)
        else:
            print("Note: No plot to show. Suggests necesesary data was not retrieved.")
    else:
        print("Note: This data structure does not support plotting.")

    return fname
