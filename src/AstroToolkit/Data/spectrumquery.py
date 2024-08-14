class SpectrumStruct(object):
    """SpectrumStruct()
    This structure is returned from spectrum queries, when read from a data file that was originally created by an spectrum query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "spectrum"

    survey: str
        survey to which the data pertains, see spectrum query description for details

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
        Returned data in the format:

        .. code-block:: python

            "wavelength": wavelength
            "flux": flux

        where:

        :**wavelength**: list<float>, wavelength of each detection
        :**flux**: flux values of each detection

    |

    """

    def __init__(self, survey, source, pos, data, identifier=None):
        self.kind = "spectrum"
        self.survey = survey
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.figure = None
        self.dataname = None

    def __str__(self):
        return "<ATK Spectrum Structure>"

    def plot(self):
        """
        Plots spectrum data and assigns the resulting figure to the **figure** attribute.

        |

        """

        from ..Plotting.plotmap import map_to_plot

        return map_to_plot(self)

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

    def showdata(self, raw=False):
        """
        Prints data structure to stdout in a readable format.

        :param bool raw: if True, doesn't collapse arrays to improve readability, default = False
        :return: self
        :rtype: class

        |

        """

        from .data_printing import print_data

        print_data(self, raw)
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
        Opens the figure stored in the **figure** attribute in the default web browser, and saves it to local files.

        :param str fname: file name to save the figure to, default file name = plotname attribute
        :return: file name to which the figure was saved
        :rtype: str

        |

        """

        from ..Plotting.plotmap import saveplot

        fname = saveplot(self, fname=fname)
        return fname


class SurveyMap(object):
    """Base class for spectrum queries"""

    def __init__(self, survey, radius, pos=None):
        self.survey = survey
        self.pos = pos
        self.radius = radius

    def query(self):
        data = globals()[f"{self.survey}_query"](pos=self.pos, radius=self.radius)
        return data


def sdss_query(pos, radius):
    from astropy import coordinates as coords
    from astropy import units as u
    from astroquery.sdss import SDSS

    position = coords.SkyCoord(pos[0], pos[1], unit="deg")
    radius = radius / 3600 * u.deg

    data = SDSS.get_spectra(coordinates=position, radius=radius, timeout=120)
    if data:
        data = data[0][1].data
        wavelength = 10 ** data["loglam"]
        flux = data["flux"] * 10**-17
        return {"wavelength": list(wavelength), "flux": list(flux)}
    else:
        print("Note: SDSS spectrum query returned no data.")
        return None


def query(survey, radius, pos=None, source=None):
    if source:
        from ..Tools import correctpm

        pos = correctpm(target_survey="sdss", source=source)
        if not pos:
            return SpectrumStruct(survey=survey, source=source, pos=pos, data=None)

    data = SurveyMap(survey=survey, radius=radius, pos=pos).query()
    return SpectrumStruct(survey=survey, source=source, pos=pos, data=data)
