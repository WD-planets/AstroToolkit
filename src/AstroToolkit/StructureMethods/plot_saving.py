from bokeh.io import output_file

newline = "\n"


def saveplot(struct, fname=None):
    from bokeh.plotting import save as bokeh_save

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()

    if not struct.figure:
        print(
            "Note: No plot was found. One will be generated with a default configuration."
        )
        struct.plot()

    if struct.figure:
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
                bokeh_save(struct.figure)
            else:
                print("Note: No plot to save.")
        else:
            print("Note: No plot to save. Create one using .plot()")
    else:
        print("Note: No plot to save. Create one using .plot()")

    return fname
