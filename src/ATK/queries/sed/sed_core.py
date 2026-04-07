import numpy as np

# AB zero point in Jy
AB_ZP_JY = 3631.0

# lambda_ref in Å, zp_vega in Jy (only for Vega-based surveys)
SED_INFO = {
    "gaia": {
        "mag_names": ["Gmag", "BPmag", "RPmag"],
        "err_names": ["e_Gmag", "e_BPmag", "e_RPmag"],
        "lambda_ref": [5822.39, 5035.75, 7619.96],
        "zp_vega": [3228.75, 3552.01, 2554.95],
        "id": "Source",
    },
    "2mass": {
        "mag_names": ["Jmag", "Hmag", "Kmag"],
        "err_names": ["e_Jmag", "e_Hmag", "e_Kmag"],
        "lambda_ref": [12350.0, 16620.0, 21590.0],
        "zp_vega": [1594.0, 1024.0, 666.8],
        "id": "2MASS",
    },
    "wise": {
        "mag_names": ["W1mag", "W2mag", "W3mag", "W4mag"],
        "err_names": ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"],
        "lambda_ref": [33526.0, 46028.0, 115608.0, 220883.0],
        "zp_vega": [309.54, 171.79, 31.67, 8.36],
        "id": "WISE",
    },
    "panstarrs": {
        "mag_names": ["gmag", "rmag", "imag", "zmag", "ymag"],
        "err_names": ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"],
        "lambda_ref": [4849.11, 6201.20, 7534.96, 8674.20, 9627.79],
        "id": "objID",
    },
    "sdss": {
        "mag_names": ["upmag", "gpmag", "rpmag", "ipmag", "zpmag"],
        "err_names": ["e_upmag", "e_gpmag", "e_rpmag", "e_ipmag", "e_zpmag"],
        "lambda_ref": [3556.52, 4702.50, 6175.58, 7489.98, 8946.71],
        "id": "objID",
    },
    "skymapper": {
        "mag_names": ["uPSF", "vPSF", "gPSF", "rPSF", "iPSF", "zPSF"],
        "err_names": ["e_uPSF", "e_vPSF", "e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF"],
        "lambda_ref": [3493.36, 3835.93, 5075.19, 6138.44, 7767.98, 9145.99],
        "id": "ObjectId",
    },
    "galex": {"mag_names": ["FUVmag", "NUVmag"], "err_names": ["e_FUVmag", "e_NUVmag"], "lambda_ref": [1535.08, 2300.78], "id": "objid"},
}


def ab_mag_to_flux_mjy(mag: float) -> float:
    """
    Convert AB magnitude to flux density in mJy
    """

    return AB_ZP_JY * 10.0 ** (-0.4 * mag) * 1e3


def get_ab_mag_offset(vega_zero_point_jy: float) -> float:
    """
    converts Vega mangitude to AB magnitude offset given Vega zero point in Jy
    """

    return -2.5 * np.log10(vega_zero_point_jy / AB_ZP_JY)
