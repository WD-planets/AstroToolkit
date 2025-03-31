class ToolInfo(object):
    def __init__(self):
        self.supported_query_kinds = ["data", "bulkdata", "reddening", "image", "lightcurve", "hrd", "sed", "spectrum"]


class SurveyInfo(object):
    @property
    def defaultSurveyTimes(self):
        data = {}
        for survey in self.defaultDataSurveys:
            data[survey] = None
        for survey in self.defaultLightcurveSurveys:
            data[survey] = None

        # data surveys
        data["gaia"] = [2016, 0]
        data["panstarrs"] = [2012, 0]
        data["skymapper"] = [2016, 0]
        data["galex"] = [2006, 8]
        data["rosat"] = [1991, 0]
        data["sdss"] = [2017, 0]
        data["wise"] = [2010, 5]
        data["twomass"] = [1999, 0]
        data["erosita"] = [2022, 0]

        # lightcurve surveys
        data["ztf"] = [2019, 0]
        data["atlas"] = [2021, 0]
        data["gaia_lc"] = [2016, 0]
        data["asassn"] = [2015, 0]
        data["crts"] = [2000, 0]
        data["tess"] = [2020, 0]

        return data

    @property
    def defaultDataSurveys(self):
        return ["gaia", "gaia_lc", "panstarrs", "skymapper", "galex", "rosat", "sdss", "wise", "twomass", "erosita"]

    @property
    def defaultReddeningSurveys(self):
        return ["stilism", "gdre"]

    @property
    def defaultLightcurveSurveys(self):
        return ["ztf", "atlas", "gaia_lc", "asassn", "crts", "tess"]

    @property
    def defaultSpectrumSurveys(self):
        return ["sdss"]

    @property
    def defaultImageSurveys(self):
        return ["panstarrs", "skymapper", "dss"]

    @property
    def dataSurveyInfo(self):
        data = {}
        for survey in self.defaultDataSurveys:
            data[survey] = {}

        data["gaia"]["mags"] = ["phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"]
        data["gaia"]["errors"] = [f"{x}_error" for x in data["gaia"]["mags"]]
        data["gaia"]["filter_wavelengths"] = [5850.88, 5041.61, 7690.74]
        data["gaia"]["catalogue"] = "I/355/gaiadr3"

        data["galex"]["mags"] = ["FUVmag", "NUVmag"]
        data["galex"]["errors"] = [f"e_{x}" for x in data["galex"]["mags"]]
        data["galex"]["filter_wavelengths"] = [2303.37, 1548.85]
        data["galex"]["catalogue"] = "II/355/galex_ais"

        data["sdss"]["mags"] = ["uPmag", "gPmag", "rPmag", "iPmag", "zPmag"]
        data["sdss"]["errors"] = [f"e_{x}" for x in data["sdss"]["mags"]]
        data["sdss"]["filter_wavelengths"] = [3608.04, 4671.78, 6141.12, 7457.89, 8922.78]
        data["sdss"]["catalogue"] = "V/154/sdss16"

        data["twomass"]["mags"] = ["Jmag", "Hmag", "Kmag"]
        data["twomass"]["errors"] = [f"e_{x}" for x in data["twomass"]["mags"]]
        data["twomass"]["filter_wavelengths"] = [12350.00, 16620.00, 21590.00]
        data["twomass"]["catalogue"] = "II/246/out"

        data["wise"]["mags"] = ["W1mag", "W2mag", "W3mag", "W4mag"]
        data["wise"]["errors"] = [f"e_{x}" for x in data["wise"]["mags"]]
        data["wise"]["filter_wavelengths"] = [33526.00, 46028.00, 115608.00, 220883.00]
        data["wise"]["catalogue"] = "II/311/wise"

        data["panstarrs"]["mags"] = ["gmag", "rmag", "imag", "zmag", "ymag"]
        data["panstarrs"]["errors"] = [f"e_{x}" for x in data["panstarrs"]["mags"]]
        data["panstarrs"]["filter_wavelengths"] = [4810.16, 6155.47, 7503.03, 8668.36, 9613.60]
        data["panstarrs"]["catalogue"] = "II/349/ps1"

        data["skymapper"]["mags"] = ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"]
        data["skymapper"]["errors"] = [f"e_{x}" for x in data["skymapper"]["mags"]]
        data["skymapper"]["filter_wavelengths"] = [5016.05, 6076.85, 6076.85, 9120.25, 3500.22, 3878.68]
        data["skymapper"]["catalogue"] = "II/379/smssdr4"

        data["rosat"]["catalogue"] = "IX/11/rosatsrc"

        data["erosita"]["catalogue"] = "J/A+A/682/A34/erass1-m"

        data["gaia_lc"]["catalogue"] = "I/355/epphot"

        return data

    @property
    def magSurveys(self):
        data = {}
        for survey, info in self.dataSurveyInfo.items():
            if "mags" in info:
                data[survey] = info

        return data

    @property
    def nonMagSurveys(self):
        data = {}
        for survey, info in self.dataSurveyInfo.items():
            if "mags" not in info:
                data[survey] = info

        return data

    @property
    def getCatalogueDict(self):
        data = self.dataSurveyInfo

        catalogueDict = {}
        for survey, info in data.items():
            catalogueDict[survey] = data[survey]["catalogue"]
        return catalogueDict

    @property
    def lightcurveSurveyInfo(self):
        data = {}
        for survey in self.defaultLightcurveSurveys:
            data[survey] = {}

        data["ztf"]["bands"] = ["g", "r", "i"]

        data["gaia_lc"]["bands"] = ["g", "bp", "rp"]

        data["asassn"]["bands"] = ["g", "v"]

        data["crts"]["bands"] = ["v"]

        data["tess"]["bands"] = ["TESS mag"]

        return data


class metadataInfo(object):
    @property
    def metadataDefaults(self):
        data = {}
        for survey in SurveyInfo().defaultDataSurveys:
            data[survey] = {}

        data["gaia"]["parameters"] = [
            "source_id",
            "ra",
            "dec",
            "pmra",
            "pmdec",
            "parallax",
            "phot_g_mean_mag",
            "phot_bp_mean_mag",
            "phot_rp_mean_mag",
        ]
        data["gaia"]["errors"] = [
            None,
            "ra_error",
            "dec_error",
            "pmra_error",
            "pmdec_error",
            "parallax_error",
            "phot_g_mean_mag_error",
            "phot_bp_mean_mag_error",
            "phot_rp_mean_mag_error",
        ]
        data["gaia"]["notes"] = [
            "source id",
            "right ascension [deg]",
            "declination [deg]",
            "Proper motion in RA [mas/yr]",
            "Proper motion in DEC [mas/yr]",
            "parallax [mas]",
            "g mag",
            "bp mag",
            "rp mag",
        ]

        data["panstarrs"]["parameters"] = ["gmag", "rmag", "imag", "zmag", "ymag"]
        data["panstarrs"]["errors"] = ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"]
        data["panstarrs"]["notes"] = ["g mag", "r mag", "i mag", "z mag", "y mag"]

        data["skymapper"]["parameters"] = ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"]
        data["skymapper"]["errors"] = ["e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF", "e_uPSF", "e_vPSF"]
        data["skymapper"]["notes"] = ["g mag", "r mag", "i mag", "z mag", "u mag", "v mag"]

        data["galex"]["parameters"] = ["NUVmag", "FUVmag"]
        data["galex"]["errors"] = ["e_NUVmag", "e_FUVmag"]
        data["galex"]["notes"] = ["FUV mag", "NUV mag"]

        data["sdss"]["parameters"] = ["gPmag", "rPmag", "iPmag", "zPmag", "uPmag"]
        data["sdss"]["errors"] = ["e_gPmag", "e_rPmag", "e_iPmag", "e_zPmag", "e_uPmag"]
        data["sdss"]["notes"] = ["g mag", "r mag", "i mag", "z mag", "u mag"]

        data["wise"]["parameters"] = ["W1mag", "W2mag", "W3mag", "W4mag"]
        data["wise"]["errors"] = ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"]
        data["wise"]["notes"] = ["W1 mag", "W2 mag", "W3 mag", "W4 mag"]

        data["twomass"]["parameters"] = ["Jmag", "Hmag", "Kmag"]
        data["twomass"]["errors"] = ["e_Jmag", "e_Hmag", "e_Kmag"]
        data["twomass"]["notes"] = ["J mag", "H mag", "K mag"]

        data["rosat"]["parameters"] = ["Name"]
        data["rosat"]["errors"] = [None]
        data["rosat"]["notes"] = ["ROSAT source name"]

        return data


class OverlayInfo(object):
    @property
    def defaultOverlayParams(self):
        data = {}
        vizierInfo = SurveyInfo().dataSurveyInfo

        # scaled detection overlays
        data["gaia"] = {"overlay_type": "detection_mag", "ra_name": "ra", "dec_name": "dec", "id_name": "designation"}
        data["galex"] = {
            "overlay_type": "detection_mag",
            "ra_name": "RAJ2000",
            "dec_name": "DEJ2000",
            "id_name": "Name",
        }
        data["wise"] = {"overlay_type": "detection_mag", "ra_name": "RAJ2000", "dec_name": "DEJ2000", "id_name": "WISE"}
        data["sdss"] = {
            "overlay_type": "detection_mag",
            "ra_name": "RA_ICRS",
            "dec_name": "DE_ICRS",
            "id_name": "objID",
        }
        data["twomass"] = {
            "overlay_type": "detection_mag",
            "ra_name": "RAJ2000",
            "dec_name": "DEJ2000",
            "id_name": "_2MASS",
        }
        data["skymapper"] = {"overlay_type": "detection_mag", "ra_name": "RAICRS", "dec_name": "DEICRS"}
        data["panstarrs"] = {
            "overlay_type": "detection_mag",
            "ra_name": "RAJ2000",
            "dec_name": "DEJ2000",
            "id_name": "objID",
        }

        # non-scaled detection overlays
        data["rosat"] = {"overlay_type": "detection", "ra_name": "RAJ2000", "dec_name": "RA_ICRS", "id_name": "Name"}
        data["erosita"] = {
            "overlay_type": "detection",
            "ra_name": "RAJ2000",
            "dec_name": "DE_ICRS",
            "id_name": "IAUName",
        }

        # tracer overlays
        data["ztf"] = {"overlay_type": "tracer"}
        data["atlas"] = {"overlay_type": "tracer"}
        data["gaia_lc"] = {"overlay_type": "tracer"}
        data["asassn"] = {"overlay_type": "tracer"}
        data["crts"] = {"overlay_type": "tracer"}

        for survey in data:
            if data[survey]["overlay_type"] == "detection_mag":
                data[survey]["mag_names"] = vizierInfo[survey]["mags"]

        return data

    @property
    def supportedOverlays(self):
        return list(self.defaultOverlayParams.keys())

    @property
    def detection_magSurveys(self):
        data = {}
        for survey, info in self.defaultOverlayParams.items():
            if info["overlay_type"] == "detection_mag":
                data[survey] = info

        return data

    @property
    def detectionSurveys(self):
        data = {}
        for survey, info in self.defaultOverlayParams.items():
            if info["overlay_type"] == "detection":
                data[survey] = info

        return data

    @property
    def tracerSurveys(self):
        data = {}
        for survey, info in self.defaultOverlayParams.items():
            if info["overlay_type"] == "tracer":
                data[survey] = info

        return data
