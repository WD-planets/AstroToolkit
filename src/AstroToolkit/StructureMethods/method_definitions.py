from ..Configuration.baseconfig import ConfigStruct
from ..Input.input_validation import check_inputs

config = ConfigStruct()
config.read_config()

newline = "\n"


def showdata(self, pprint, print_methods):
    """
    Prints data structure to stdout in a readable format.

    :param pprint: collapse arrays and other objects to improve readability, defaults to True
    :type pprint: bool, optional

    :return: Self

    |

    """

    corrected_inputs = check_inputs(
        {"pprint": [pprint, bool], "print_methods": [print_methods, bool]}, label="showdata"
    )
    pprint, print_methods = corrected_inputs

    from .data_printing import print_data

    print_data(self, pprint, print_methods)

    return self


def savedata(self, fname):
    """
    Saves a data structure's data to local files.

    :param fname: overrides file name, defaults to file name given by the data structure's 'dataname' attribute
    :type fname: str, optional

    :return: name of file to which data was saved
    :rtype: str

    |

    """
    corrected_inputs = check_inputs({"fname": [fname, str]}, label="savedata")
    fname = corrected_inputs[0]

    from .data_saving import savedata

    fname = savedata(self, fname)

    return fname


def showplot(self, fname):
    """
    Opens the figure stored in the 'figure' attribute in the default web browser, and saves it to local files. If the data structure doesn't yet hold a figure, one will be generated with a default configuration first.

    :param fname: file name to save the figure to, defaults to file name given by the data structure's 'plotname' attribute
    :type fname: str, optional

    :return: file name to which the figure was saved
    :rtype: str

    |

    """

    from .plot_showing import showplot

    corrected_inputs = check_inputs({"fname": [fname, str]}, label="showplot")
    fname = corrected_inputs[0]

    fname = showplot(self, fname)
    return fname


def saveplot(self, fname):
    """
    Saves the figure stored in the 'figure' attribute to local files without opening it in the web browser. If the data structure doesn't yet hold a figure, one will be generated with a default configuration first.

    :param fname: file name to save the figure to, defaults to file name given by the data structure's 'plotname' attribute
    :type fname: str, optional

    :return: file name to which the figure was saved
    :rtype: str

    |

    """

    corrected_inputs = check_inputs({"fname": [fname, str]}, label="saveplot")
    fname = corrected_inputs[0]

    from .plot_saving import saveplot

    fname = saveplot(self, fname=fname)
    return fname


def exportplot(self, fname):
    """
    Exports the figure stored in the 'figure' attribute to a PNG, and saves it to local files. If the data structure doesn't yet hold a figure, one will be generated with a default configuration first.

    :param fname: file name to same the figure to, defaults to file name given by the data structure's 'plotname' attribute
    :type fname: str, optional

    :return: file name to which the figure was saved
    :rtype: str

    |

    """

    corrected_inputs = check_inputs({"fname": [fname, str]}, "exportplot")
    fname = corrected_inputs[0]

    from .plot_exporting import exportplot

    fname = exportplot(self, fname)
    return fname


def plot(self, **kwargs):
    """
    Plots data contained within a given data stucture and assigns the resulting figure to the data structure's 'figure' attribute.

    **Note:** Some data structures may allow for additional optional arguments, see their specific plot() methods for details.

    |

    """

    if config.enable_notifications:
        print(f"Plotting {self.kind} data...{newline}")

    from ..Plotting.plotting_map import get_plot

    defaults = {
        "kind": ["lightcurve", str],
        "colours": [None, list],
        "bands": [None, list],
        "spectrum_overlay": [None, None],
        "survey": ["sdss", str],
        "freq": [None, float],
        "bins": [None, int],
        "timeformat": ["reduced", str],
        "method": ["ls", str],
        "foverlay": [True, bool],
        "repeat": [2, int],
        "shift": [0, float],
        "start_freq": [0, float],
        "stop_freq": [60, float],
        "samples": [150000, int],
        "searchradius": [float(config.overlay_simbad_search_radius), float],
    }

    if self.kind == "lightcurve":
        if "kind" in kwargs:
            plot_kind = kwargs["kind"]
        else:
            plot_kind = "lightcurve"

        if plot_kind == "lightcurve":
            parameters = ["kind", "colours", "bands", "timeformat"]
        elif plot_kind == "powspec":
            parameters = ["kind", "method", "start_freq", "stop_freq", "samples"]
        elif plot_kind == "phasefold":
            parameters = ["kind", "freq", "bins", "method", "foverlay", "repeat", "shift"]
    elif self.kind == "image":
        parameters = ["searchradius"]
    elif self.kind == "sed":
        parameters = ["spectrum_overlay"]
    else:
        parameters = []

    inputs = {}
    for parameter in parameters:
        if parameter not in kwargs:
            kwargs[parameter] = defaults[parameter][0]
        inputs[parameter] = [kwargs[parameter], defaults[parameter][1]]
    inputs["data_kind"] = [self.kind, str]

    corrected_inputs = check_inputs(inputs, label="plot")
    for index, parameter in enumerate(parameters):
        kwargs[parameter] = corrected_inputs[index]

    return get_plot(self, **kwargs)
