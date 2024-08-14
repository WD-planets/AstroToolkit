import warnings

import astropy.coordinates as coord
import astropy.units as u
import pandas as pd
from astroquery.vizier import Vizier

from ..Misc.pmcorrection import correctpm
from .rename_headers import renameHeadersDR3

warnings.simplefilter(action="ignore", category=UserWarning)

# ensure that the row limit in returned data is infinite
row_limit = -1
Vizier.ROW_LIMIT = -1


class DataStruct(object):
    """DataStruct()
    This structure is returned from data, phot, bulkphot and reddening queries, when read from a data file that was originally created by a data, phot, bulkphot or reddening query, or through the Models module (in which case all attributes are set to None).

    Attributes
    ----------

    kind : str
        "data"

    subkind: str
        query kind, from: "data", "phot", "bulkphot", "reddening"

    survey: str
        survey to which the data pertains, see data query description for details. For bulkphot queries, survey = None.

    catalogue: str
        Vizier catalogue ID

    source: int
        Gaia source ID (if one was used to target the system, else None)

    pos: list<int>
        Target position [right ascension, declination] in degrees

    identifier: str
        Target position in JHHMMSS.SS±DDMMSS.SS format

    dataname: str
        Default file name that data will be saved to using the savedata() method

    data: dict<float>
        Returned data with keys = column headers, items = row values

    |

    """

    def __init__(
        self,
        survey,
        catalogue,
        source,
        pos,
        data,
        identifier=None,
        sub_kind="data",
    ):
        self.kind = "data"
        self.subkind = sub_kind
        self.survey = survey
        self.catalogue = catalogue
        self.source = source
        self.pos = pos
        self.identifier = identifier
        self.data = data
        self.dataname = None

    def __str__(self):
        return "<ATK Data Structure>"

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

        :param bool raw: collapse arrays to improve readability, default = True
        :return: self
        :rtype: class

        |

        """

        from .data_printing import print_data

        print_data(self, raw)
        return self


class SurveyInfo:
    def __init__(self):
        self.list = [
            "gaia",
            "gaia_lc",
            "panstarrs",
            "skymapper",
            "galex",
            "rosat",
            "sdss",
            "wise",
            "twomass",
            "erosita",
        ]
        self.times = {
            "gaia": [2016, 0],
            "panstarrs": [2012, 0],
            "skymapper": [2016, 0],
            "galex": [2006, 8],
            "rosat": [1991, 0],
            "sdss": [2017, 0],
            "wise": [2010, 5],
            "twomass": [1999, 0],
            "ztf": [2019, 0],
            "erosita": [2022, 0],
            "atlas": [2021, 0],
            "gaia_lc": [2016, 0],
            "asassn": [2015, 0],
            "crts": [2000, 0],
            "tess": [2020, 0],
        }
        self.catalogues = {
            "gaia": "I/355/gaiadr3",
            "panstarrs": "II/349/ps1",
            "skymapper": "II/379/smssdr4",
            "galex": "II/335/galex_ais",
            "rosat": "IX/11/rosatsrc",
            "sdss": "V/154/sdss16",
            "wise": "II/311/wise",
            "twomass": "II/246/out",
            "erosita": "J/A+A/682/A34/erass1-m",
            "gaia_lc": "I/355/epphot",
        }

        self.bulkphot_surveys = [
            "gaia",
            "galex",
            "sdss",
            "twomass",
            "wise",
            "panstarrs",
            "skymapper",
        ]

        self.sed_param_names = {
            "gaia": {
                "filter_wavelengths": [5850.88, 5041.61, 7690.74],
                "mag_names": [
                    "phot_g_mean_mag",
                    "phot_bp_mean_mag",
                    "phot_rp_mean_mag",
                ],
                "error_names": [
                    "phot_g_mean_mag_error",
                    "phot_bp_mean_mag_error",
                    "phot_rp_mean_mag_error",
                ],
            },
            "galex": {
                "filter_wavelengths": [2303.37, 1548.85],
                "mag_names": ["NUVmag", "FUVmag"],
                "error_names": ["e_NUVmag", "e_FUVmag"],
            },
            "sdss": {
                "filter_wavelengths": [3608.04, 4671.78, 6141.12, 7457.89, 8922.78],
                "mag_names": ["uPmag", "gPmag", "rPmag", "iPmag", "zPmag"],
                "error_names": ["e_uPmag", "e_gPmag", "e_rPmag", "e_iPmag", "e_zPmag"],
            },
            "twomass": {
                "filter_wavelengths": [12350.00, 16620.00, 21590.00],
                "mag_names": ["Jmag", "Hmag", "Kmag"],
                "error_names": ["e_Jmag", "e_Hmag", "e_Kmag"],
            },
            "wise": {
                "filter_wavelengths": [33526.00, 46028.00, 115608.00, 220883.00],
                "mag_names": ["W1mag", "W2mag", "W3mag", "W4mag"],
                "error_names": ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"],
            },
            "panstarrs": {
                "filter_wavelengths": [4810.16, 6155.47, 7503.03, 8668.36, 9613.60],
                "mag_names": ["gmag", "rmag", "imag", "zmag", "ymag"],
                "error_names": ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"],
            },
            "skymapper": {
                "filter_wavelengths": [
                    5016.05,
                    6076.85,
                    6076.85,
                    9120.25,
                    3500.22,
                    3878.68,
                ],
                "mag_names": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
                "error_names": [
                    "e_gPSF",
                    "e_rPSF",
                    "e_iPSF",
                    "e_zPSF",
                    "e_uPSF",
                    "e_vPSF",
                ],
            },
        }

        self.lightcurve_bands = {
            "ztf": [
                "g",
                "r",
                "i",
            ],
            "atlas": ["o", "c", "i"],
            "gaia": ["g", "bp", "rp"],
            "asassn": ["g", "v"],
            "crts": ["v"],
            "tess": ["i"],
        }

        self.metadata_defaults = {
            "gaia": {
                "parameters": [
                    "source_id",
                    "ra",
                    "dec",
                    "pmra",
                    "pmdec",
                    "parallax",
                    "phot_g_mean_mag",
                    "phot_bp_mean_mag",
                    "phot_rp_mean_mag",
                ],
                "errors": [
                    None,
                    "ra_error",
                    "dec_error",
                    "pmra_error",
                    "pmdec_error",
                    "parallax_error",
                    "phot_g_mean_mag_error",
                    "phot_bp_mean_mag_error",
                    "phot_rp_mean_mag_error",
                ],
                "notes": [
                    "source id",
                    "right ascension [deg]",
                    "declination [deg]",
                    "Proper motion in RA direction [mas/yr]",
                    "Proper motion in DEC direction [mas/yr]",
                    "parallax [mas]",
                    "g mag",
                    "bp mag",
                    "rp mag",
                ],
            },
            "panstarrs": {
                "parameters": [
                    "gmag",
                    "rmag",
                    "imag",
                    "zmag",
                    "ymag",
                ],
                "errors": [
                    "e_gmag",
                    "e_rmag",
                    "e_imag",
                    "e_zmag",
                    "e_ymag",
                ],
                "notes": ["g mag", "r mag", "i mag", "z mag", "y mag"],
            },
            "skymapper": {
                "parameters": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
                "errors": [
                    "e_gPSF",
                    "e_rPSF",
                    "e_iPSF",
                    "e_zPSF",
                    "e_uPSF",
                    "e_vPSF",
                ],
                "notes": ["g mag", "r mag", "i mag", "z mag", "u mag", "v mag"],
            },
            "galex": {
                "parameters": ["NUVmag", "FUVmag"],
                "errors": ["e_NUVmag", "e_FUVmag"],
                "notes": ["FUV mag", "NUV mag"],
            },
            "sdss": {
                "parameters": ["gPmag", "rPmag", "iPmag", "zPmag", "uPmag"],
                "errors": ["e_gPmag", "e_rPmag", "e_iPmag", "e_zPmag", "e_uPmag"],
                "notes": ["g mag", "r mag", "i mag", "z mag", "u mag"],
            },
            "wise": {
                "parameters": ["W1mag", "W2mag", "W3mag", "W4mag"],
                "errors": ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"],
                "notes": ["W1 mag", "W2 mag", "W3 mag", "W4 mag"],
            },
            "twomass": {
                "parameters": ["Jmag", "Hmag", "Kmag"],
                "errors": ["e_Jmag", "e_Hmag", "e_Kmag"],
                "notes": ["J mag", "H mag", "K mag"],
            },
            "rosat": {
                "parameters": ["Name"],
                "errors": [None],
                "notes": ["ROSAT source name"],
            },
        }

        self.spectrum_surveys = ["sdss"]

        self.image_surveys = ["panstarrs", "skymapper", "dss", "any"]

        self.supported_overlays = [
            "gaia",
            "galex",
            "wise",
            "sdss",
            "twomass",
            "skymapper",
            "panstarrs",
            "rosat",
            "erosita",
            "atlas",
            "gaia_lc",
            "asassn",
            "crts",
            "ztf",
        ]

        self.lightcurve_surveys = ["atlas", "ztf", "crts", "asassn", "gaia", "tess"]

        self.reddening_surveys = ["stilism", "gdre"]

        self.supported_query_kinds = [
            "data",
            "phot",
            "bulkphot",
            "reddening",
            "image",
            "lightcurve",
            "hrd",
            "sed",
            "spectrum",
        ]

    @property
    def overlay_param_names(self):
        data = {
            "gaia": {
                "overlay_type": "detection_mag",
                "ra_name": "ra",
                "dec_name": "dec",
                "colours": ["limegreen", "blue", "red"],
                "default_overlay_mag": "phot_g_mean_mag",
            },
            "galex": {
                "overlay_type": "detection_mag",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colours": ["purple", "violet"],
                "default_overlay_mag": "NUVmag",
            },
            "wise": {
                "overlay_type": "detection_mag",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colours": ["firebrick", "orange", "gold", "yellow"],
                "default_overlay_mag": "W1mag",
            },
            "sdss": {
                "overlay_type": "detection_mag",
                "ra_name": "RA_ICRS",
                "dec_name": "DE_ICRS",
                "colours": ["tomato", "darkorange", "khaki", "aqua", "mediumblue"],
                "default_overlay_mag": "gPmag",
            },
            "twomass": {
                "overlay_type": "detection_mag",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colours": ["orangered", "goldenrod", "lightyellow"],
                "default_overlay_mag": "Jmag",
            },
            "skymapper": {
                "overlay_type": "detection_mag",
                "ra_name": "RAICRS",
                "dec_name": "DEICRS",
                "colours": [
                    "indianred",
                    "darkgoldenrod",
                    "lawngreen",
                    "dodgerblue",
                    "stateblue",
                    "blueviolet",
                ],
                "default_overlay_mag": "gPSF",
            },
            "panstarrs": {
                "overlay_type": "detection_mag",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colours": [
                    "salmon",
                    "yellowgreen",
                    "turquoise",
                    "midnightblue",
                    "indigo",
                ],
                "default_overlay_mag": "gmag",
            },
            "rosat": {
                "overlay_type": "detection",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colour": "deeppink",
            },
            "erosita": {
                "overlay_type": "detection",
                "ra_name": "RA_ICRS",
                "dec_name": "DE_ICRS",
                "colour": "lightpink",
            },
            "atlas": {"overlay_type": "tracer", "colour": "bisque"},
            "gaia_lc": {"overlay_type": "tracer", "colour": "forestgreen"},
            "asassn": {"overlay_type": "tracer", "colour": "mediumorchid"},
            "crts": {"overlay_type": "tracer", "colour": "teal"},
            "ztf": {"overlay_type": "tracer", "colour": "lightsteelblue"},
        }

        for survey in data:
            if data[survey]["overlay_type"] == "detection_mag":
                data[survey]["mag_names"] = self.sed_param_names[survey]["mag_names"]
                data[survey]["marker_type"] = "circle"
            elif (
                data[survey]["overlay_type"] == "detection"
                or data[survey]["overlay_type"] == "tracer"
            ):
                data[survey]["marker_type"] = "cross"

        return data


class VizierQuery(object):
    """Performs Vizier queries"""

    def __init__(self, catalogue, radius, survey=None, pos=None, source=None):
        self.catalogue = catalogue
        self.pos = pos
        self.radius = radius
        self.survey = survey
        self.source = source
        self.data = None

        self.f_return = DataStruct(
            survey=self.survey,
            catalogue=self.catalogue,
            source=self.source,
            pos=self.pos,
            data=None,
        )

    # performs queries by coordinates
    def pos_query(self):
        data = []
        v = Vizier(columns=["**"], row_limit=row_limit)
        data.append(
            v.query_region(
                coord.SkyCoord(
                    ra=self.pos[0], dec=self.pos[1], unit=(u.deg, u.deg), frame="icrs"
                ),
                width=self.radius * u.arcsec,
                catalog=self.catalogue,
            )
        )
        try:
            data = data[0][0].to_pandas().reset_index(drop=True)
            return self.check_data(data)
        except:
            return self.f_return

    # performs queries by Gaia source
    def source_query(self):
        v = Vizier(
            columns=["**"],
            column_filters={"Source": "==" + str(self.source)},
            row_limit=row_limit,
        )
        data = v.query_object(f"GAIA DR3 {self.source}", catalog=self.catalogue)
        try:
            data = data[0].to_pandas().reset_index(drop=True)
            return self.check_data(data)
        except:
            return self.f_return

    # checks if any data was returned
    def check_data(self, data):
        try:
            if not data.empty:
                data = pd.DataFrame.to_dict(data, orient="list")
                return DataStruct(
                    survey=self.survey,
                    catalogue=self.catalogue,
                    source=self.source,
                    pos=self.pos,
                    data=data,
                )
            else:
                return DataStruct(
                    survey=self.survey,
                    catalogue=self.catalogue,
                    source=self.source,
                    pos=self.pos,
                    data=None,
                )
        except:
            return DataStruct(
                survey=self.survey,
                catalogue=self.catalogue,
                source=self.source,
                pos=self.pos,
                data=None,
            )


# main
# maps coordinates to vizier surveys, performing proper motion correction for source queries.
def query(survey, radius, pos=None, source=None):
    # get the necessary basic survey info
    supported_surveys = SurveyInfo().list
    supported_catalogues = SurveyInfo().catalogues
    survey_times = SurveyInfo().times

    # if survey isn't a supported survey, take the 'survey' to be a vizier catalogue ID
    if survey not in supported_surveys:
        catalogue = survey
    else:
        catalogue = supported_catalogues[survey]

    # perform coordinate Vizier query
    if pos:
        data = VizierQuery(catalogue=catalogue, radius=radius, pos=pos).pos_query()

    # perform source Vizier query
    elif source:
        if catalogue == "I/355/gaiadr3":
            data = VizierQuery(
                survey=survey, catalogue=catalogue, radius=radius, source=source
            ).source_query()
        elif catalogue == "I/355/epphot":
            data = VizierQuery(
                survey=survey, catalogue=catalogue, radius=radius, source=source
            ).source_query()
        else:
            gaia_data = (
                VizierQuery(
                    survey=survey,
                    catalogue=supported_catalogues["gaia"],
                    radius=radius,
                    source=source,
                )
                .source_query()
                .data
            )
            if gaia_data:
                ra, dec, pmra, pmdec = (
                    gaia_data["RA_ICRS"][0],
                    gaia_data["DE_ICRS"][0],
                    gaia_data["pmRA"][0],
                    gaia_data["pmDE"][0],
                )
            else:
                print("Note: data query unsuccessful or returned no data.")
                return DataStruct(
                    survey=survey,
                    catalogue=catalogue,
                    source=source,
                    pos=pos,
                    data=None,
                )

            if survey in supported_surveys:
                ra, dec = correctpm(
                    [2016, 0], survey_times[str(survey)], ra, dec, pmra, pmdec
                )
            data = VizierQuery(
                survey=survey,
                catalogue=catalogue,
                radius=radius,
                pos=[ra, dec],
                source=source,
            ).pos_query()
    else:
        raise Exception("No source or coordinates provided.")

    if catalogue == "I/355/gaiadr3" and data.data:
        data = renameHeadersDR3(data)

    if not data.data:
        print(f"Note: {survey} data query unsuccessful or returned no data.")

    return data
