class SedStruct(object):
    """SedStruct()
    This structure is returned from sed queries, when read from a data file that was originally created by an sed query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "sed"

    source: int
        Gaia source ID (if one was used to target the system, else None)

    pos: list<int>
        Target position [right ascension, declination] in degrees

    identifier: str
        Target position in JHHMMSS.SS±DDMMSS.SS format

    dataname: str
        Default file name that data will be saved to when using the savedata() method

    plotname: str
        Default file name that plot will be saved to when using the showplot() and saveplot() methods

    figure: None
        Stores any figure resulting from the plot() method

    data: dict<list>
        Returned data as a dict of entries, with keys = survey, with each entry in the format:

        .. code-block:: python

            "wavelength": wavelength
            "flux": flux
            "flux_rel_err" flux_rel_err

        where:

        :**wavelength**: list<float>, central filter wavelengths of each detection
        :**flux**: flux values of each detection
        :**flux_rel_err**: relative error on the flux

    |

    """

    def __init__(self, source, pos, data, identifier=None):
        self.kind = "sed"
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.figure = None
        self.dataname = None

    def __str__(self):
        return "<ATK SED Structure>"

    def plot(self, spectrum_overlay=False, survey=None):
        """
        Plots sed data and assigns the resulting figure to the **figure** attribute.

        :param bool spectrum_overlay: if True, overlays a spectrum from chosen survey if available
        :param str survey: survey to fetch spectrum from. See spectrum query for available surveys.

        |

        """

        from ..Plotting.plotmap import map_to_plot

        return map_to_plot(self, spectrum_overlay=spectrum_overlay, survey=survey)

    def savedata(self, fname=None):
        """
        Saves data to local files. Default file name = dataname attribute.

        :param str fname: overrides file name
        :return: name of file to which data was saved
        :rtype: str

        |

        """

        from .data_saving import savedata

        fname = savedata(self, fname)
        return fname

    def showdata(self):
        """
        Prints data structure to stdout in a readable format.

        :param bool raw: if True, doesn't collapse arrays to improve readability, default = False
        :return: self
        :rtype: class

        |

        """

        from .data_printing import print_data

        print_data(self)
        return self

    def showplot(self, fname=None):
        """
        Opens the figure stored in the **figure** attribute in the default web browser, and saves it to local files.

        :param str fname: file name to save the figure to, default file name = plotname attribute
        :return: file name to which the figure was saved
        :rtype: str

        |

        """

        from ..Plotting.plotmap import showplot

        fname = showplot(self, fname=fname)
        return fname

    def saveplot(self, fname=None):
        """
        Saves the figure stored in the **figure** attribute to local files without opening it in the web browser.

        :param str fname: file name to save the figure to, default file name = plotname attribute
        :return: file name to which the figure was saved
        :rtype: str

        |

        """

        from ..Plotting.plotmap import saveplot

        fname = saveplot(self, fname=fname)
        return fname


def mag_to_flux(mag, zp, wl):
    # sometimes overflows if a bad mag is passed (e.g. -999 for some surveys)
    try:
        flux = zp * 10 ** (-0.4 * mag)
    except:
        return None
    # convert to mjy
    c = 2.988 * 10**18
    fnl = 1 * 10 ** (-23)
    flux = flux / ((fnl * c) / wl**2) * 1000

    return flux


def format_data(survey, photometry, filter_wavelengths, mag_names, error_names):
    import numpy as np

    if survey != "gaia":
        zero_points = [
            10 ** ((5 * np.log10(x) + 2.406) / -2.5) for x in filter_wavelengths
        ]
    else:
        zero_points = [2.5e-9, 4.11e-9, 1.24e-9]

    sed_datapoints = {
        "survey": survey,
        "wavelength": [],
        "flux": [],
        "flux_rel_err": [],
    }
    for filter_wavelength, mag_name, error_name, zero_point in zip(
        filter_wavelengths, mag_names, error_names, zero_points
    ):
        mag, mag_err = photometry[mag_name][0], photometry[error_name][0]

        flux = mag_to_flux(mag=mag, zp=zero_point, wl=filter_wavelength)
        if flux:
            rel_err = flux * mag_err / mag

            sed_datapoints["wavelength"].append(filter_wavelength)
            sed_datapoints["flux"].append(flux)
            sed_datapoints["flux_rel_err"].append(rel_err)

    return sed_datapoints


def query(radius, pos=None, source=None):
    from ..Data.dataquery import SurveyInfo
    from ..Tools import query

    sed_params = SurveyInfo().sed_param_names

    bulkphot = query(
        kind="bulkphot", pos=pos, source=source, radius=radius, level="internal"
    )
    if bulkphot.data:
        pos = bulkphot.data["gaia"]["ra"][0], bulkphot.data["gaia"]["dec"][0]

        bulkphot.data = {
            key: value for key, value in bulkphot.data.items() if value is not None
        }

    sed_data = []
    for survey in bulkphot.data:
        filter_wavelengths, mag_names, error_names = (
            sed_params[survey]["filter_wavelengths"],
            sed_params[survey]["mag_names"],
            sed_params[survey]["error_names"],
        )

        sed_data.append(
            format_data(
                survey=survey,
                photometry=bulkphot.data[survey],
                filter_wavelengths=filter_wavelengths,
                mag_names=mag_names,
                error_names=error_names,
            )
        )

    data_struct = SedStruct(pos=pos, source=source, data=sed_data)
    return data_struct
