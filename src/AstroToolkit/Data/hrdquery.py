class HrdStruct(object):
    """HrdStruct()
    This structure is returned from hrd queries, when read from a hrd file that was originally created by an image query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "hrd"

    survey: str
        "gaia"

    sources: list<int>
        Gaia source IDs

    identifiers: list<str>
        Target positions in JHHMMSS.SS±DDMMSS.SS format

    dataname: str
        Default file name that data will be saved to when using the savedata() method

    plotname: str
        Default file name that plot will be saved to when using the showplot() and saveplot() methods

    figure: None
        Stores any figure resulting from the plot() method

    data: dict<list>
        Returned data in format:

        .. code-block:: python

            "bp-rp": bp-rp
            "absg": absg

        where:

        :**bp-rp**: list<float>, bp-rp colour for each source in sources
        :**absg**: list<float>, absolute g magnitude for each source in sources

    |

    """

    def __init__(self, sources, data, identifiers=None, survey="gaia"):
        self.kind = "hrd"

        self.survey = survey
        self.sources = sources
        self.identifiers = identifiers
        self.data = data
        self.figure = None
        self.dataname = None

    def __str__(self):
        return "<ATK HRD Structure>"

    def plot(self):
        """
        Plots hrd data and assigns the resulting figure to the **figure** attribute.

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
        Saves data to local files. Default file name = dataname attribute.

        :param str fname: overrides file name
        :return: name of file to which data was saved
        :rtype: str

        |

        """

        from .data_printing import print_data

        print_data(self, raw)
        return self

    def showplot(self, fname=None):
        """
        Saves data to local files. Default file name = dataname attribute.

        :param str fname: overrides file name
        :return: name of file to which data was saved
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


def gather_data(sources):
    import numpy as np

    if not isinstance(sources, list):
        sources = [sources]

    x, y = [], []

    from ..Tools import query

    bad_indices = []
    for index, source in enumerate(sources):
        gaia_data = query(
            kind="data", source=source, survey="gaia", level="internal"
        ).data
        if gaia_data:
            gmag, bpmag, rpmag, parallax = (
                gaia_data["phot_g_mean_mag"][0],
                gaia_data["phot_bp_mean_mag"][0],
                gaia_data["phot_rp_mean_mag"][0],
                gaia_data["parallax"][0],
            )

            x.append(bpmag - rpmag)
            y.append(gmag + 5 * np.log10(parallax / 1000) + 5)
        else:
            bad_indices.append(index)

    sources_formatted = [
        source for i, source in enumerate(sources) if i not in bad_indices
    ]

    return HrdStruct(sources=sources_formatted, data={"bp-rp": x, "absg": y})
