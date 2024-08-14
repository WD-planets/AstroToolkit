from ..Data.dataquery import SurveyInfo

bulkphot_surveys = SurveyInfo().bulkphot_surveys


def query(survey, radius, pos=None, source=None):
    from ..Tools import query

    dataStruct = query(
        kind="data",
        survey=survey,
        pos=pos,
        source=source,
        radius=radius,
        level="internal",
    )

    dataStruct.subkind = "phot"

    if dataStruct.data:
        if survey == "gaia":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "ra",
                    "dec",
                    "source_id",
                    "phot_g_mean_mag",
                    "phot_g_mean_mag_error",
                    "phot_bp_mean_mag",
                    "phot_bp_mean_mag_error",
                    "phot_rp_mean_mag",
                    "phot_rp_mean_mag_error",
                ]
            }
        elif survey == "galex":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RAJ2000",
                    "DEJ2000",
                    "objid",
                    "NUVmag",
                    "e_NUVmag",
                    "FUVmag",
                    "e_FUVmag",
                ]
            }
        elif survey == "sdss":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RA_ICRS",
                    "DE_ICRS",
                    "objID",
                    "uPmag",
                    "e_uPmag",
                    "gPmag",
                    "e_gPmag",
                    "rPmag",
                    "e_rPmag",
                    "iPmag",
                    "e_iPmag",
                    "zPmag",
                    "e_zPmag",
                ]
            }
        elif survey == "twomass":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RAJ2000",
                    "DEJ2000",
                    "_2MASS",
                    "Jmag",
                    "e_Jmag",
                    "Hmag",
                    "e_Hmag",
                    "Kmag",
                    "e_Kmag",
                ]
            }
        elif survey == "wise":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RAJ2000",
                    "DEJ2000",
                    "WISE",
                    "W1mag",
                    "e_W1mag",
                    "W2mag",
                    "e_W2mag",
                    "W3mag",
                    "e_W3mag",
                    "W4mag",
                    "e_W4mag",
                ]
            }
        elif survey == "panstarrs":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RAJ2000",
                    "DEJ2000",
                    "objID",
                    "gmag",
                    "e_gmag",
                    "rmag",
                    "e_rmag",
                    "imag",
                    "e_imag",
                    "zmag",
                    "e_zmag",
                    "ymag",
                    "e_ymag",
                ]
            }
        elif survey == "skymapper":
            photometry = {
                key: dataStruct.data[key]
                for key in [
                    "RAICRS",
                    "DEICRS",
                    "ObjectId",
                    "gPSF",
                    "e_gPSF",
                    "rPSF",
                    "e_rPSF",
                    "iPSF",
                    "e_iPSF",
                    "zPSF",
                    "e_zPSF",
                    "uPSF",
                    "e_uPSF",
                    "vPSF",
                    "e_vPSF",
                ]
            }
        else:
            raise Exception(f"Unsupported survey. Accepted surveys: {bulkphot_surveys}")

        dataStruct.data = photometry

    return dataStruct


def bulkphot_query(radius, pos=None, source=None):
    from ..Tools import query

    bulk_phot = {}
    for survey in bulkphot_surveys:
        data = query(
            kind="phot",
            pos=pos,
            source=source,
            radius=radius,
            survey=survey,
            level="internal",
        ).data
        bulk_phot[survey] = data

    if source:
        gaia_data = query(
            kind="data", source=source, survey="gaia", level="internal"
        ).data
        if gaia_data:
            pos = [gaia_data["ra"][0], gaia_data["dec"][0]]

    from ..Data.dataquery import DataStruct

    dataStruct = DataStruct(
        survey="all",
        catalogue=None,
        pos=pos,
        source=source,
        data=bulk_phot,
        sub_kind="bulkphot",
    )

    return dataStruct
