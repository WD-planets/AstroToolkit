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

        self.sed_param_names = {
            "gaia": {
                "filter_wavelengths": [5850.88, 5041.61, 7690.74],
                "mag_names": ["phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"],
                "error_names": ["phot_g_mean_mag_error", "phot_bp_mean_mag_error", "phot_rp_mean_mag_error"],
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
                "filter_wavelengths": [5016.05, 6076.85, 6076.85, 9120.25, 3500.22, 3878.68],
                "mag_names": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
                "error_names": ["e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF", "e_uPSF", "e_vPSF"],
            },
        }

        self.lightcurve_bands = {
            "ztf": ["g", "r", "i"],
            "atlas": ["o", "c", "i"],
            "gaia": ["g", "bp", "rp"],
            "asassn": ["g", "v"],
            "crts": ["v"],
            "tess": ["TESS mag"],
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
                    "Proper motion in RA [mas/yr]",
                    "Proper motion in DEC [mas/yr]",
                    "parallax [mas]",
                    "g mag",
                    "bp mag",
                    "rp mag",
                ],
            },
            "panstarrs": {
                "parameters": ["gmag", "rmag", "imag", "zmag", "ymag"],
                "errors": ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"],
                "notes": ["g mag", "r mag", "i mag", "z mag", "y mag"],
            },
            "skymapper": {
                "parameters": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
                "errors": ["e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF", "e_uPSF", "e_vPSF"],
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
            "rosat": {"parameters": ["Name"], "errors": [None], "notes": ["ROSAT source name"]},
        }

        self.spectrum_surveys = ["sdss"]

        self.image_surveys = ["panstarrs", "skymapper", "dss"]

        self.marker_overlays = ["gaia", "galex", "wise", "sdss", "twomass", "skymapper", "panstarrs"]

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

        self.supported_query_kinds = ["data", "bulkdata", "reddening", "image", "lightcurve", "hrd", "sed", "spectrum"]

        self.survey_id_names = {
            "gaia": "designation",
            "galex": "Name",
            "panstarrs": "objID",
            "skymapper": "ObjectId",
            "rosat": "Name",
            "sdss": "objID",
            "wise": "WISE",
            "twomass": "_2MASS",
            "erosita": "IAUName",
        }

        self.time_units = {"ztf": "hjd", "crts": "mjd", "asassn": "mjd", "gaia": "mjd", "atlas": "mjd", "tess": "mjd"}

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
                "colours": ["indianred", "darkgoldenrod", "lawngreen", "dodgerblue", "stateblue", "blueviolet"],
                "default_overlay_mag": "gPSF",
            },
            "panstarrs": {
                "overlay_type": "detection_mag",
                "ra_name": "RAJ2000",
                "dec_name": "DEJ2000",
                "colours": ["salmon", "yellowgreen", "turquoise", "midnightblue", "indigo"],
                "default_overlay_mag": "gmag",
            },
            "rosat": {"overlay_type": "detection", "ra_name": "RAJ2000", "dec_name": "DEJ2000", "colour": "deeppink"},
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
            "ztf": {"overlay_type": "tracer", "colour": "orange"},
        }

        for survey in data:
            if data[survey]["overlay_type"] == "detection_mag":
                data[survey]["mag_names"] = self.sed_param_names[survey]["mag_names"]
                data[survey]["marker_type"] = "circle"
            elif data[survey]["overlay_type"] == "detection" or data[survey]["overlay_type"] == "tracer":
                data[survey]["marker_type"] = "cross"

        return data


class OverlayInfo(object):
    def __init__(self):
        self.sections = ["scaled_detection_surveys", "detection_surveys", "tracer_surveys"]
        self.scaled_detection_surveys = ["gaia", "galex", "wise", "sdss", "twomass", "skymapper", "panstarrs"]
        self.detection_surveys = ["rosat", "erosita"]
        self.tracer_surveys = ["ztf", "atlas", "crts", "asassn", "gaia_lc"]

    @property
    def defaultOverlayParams(self):
        data = {}
        for section in self.sections:
            data[section] = {}
            for survey in getattr(self, section):
                data[section][survey] = {}

        data["scaled_detection_surveys"]["gaia"] = {
            "ra_name": "ra",
            "dec_name": "dec",
            "mag_names": ["phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"],
            "default_overlay_mag": "phot_g_mean_mag",
        }
        data["scaled_detection_surveys"]["galex"] = {
            "overlay_type": "detection_mag",
            "ra_name": "RAJ2000",
            "dec_name": "DEJ2000",
            "colours": ["purple", "violet"],
            "default_overlay_mag": "NUVmag",
        }
