from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_overlays.ini"

DEFAULTS = {
    "photometric": {
        "gaia": {
            "mags": ["Gmag", "BPmag", "RPmag"],
            "errors": ["e_Gmag", "e_BPmag", "e_RPmag"],
            "lon_column": "RA_ICRS",
            "lat_column": "DE_ICRS",
            "frame": "icrs",
            "id_column": "Source",
        },
        "galex": {
            "mags": ["NUVmag", "FUVmag"],
            "errors": ["e_NUVmag", "e_FUVmag"],
            "lon_column": "RAJ2000",
            "lat_column": "DEJ2000",
            "frame": "icrs",
            "id_column": "Name",
        },
        "wise": {
            "mags": ["W1mag", "W2mag", "W3mag", "W4mag"],
            "errors": ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"],
            "lon_column": "RAJ2000",
            "lat_column": "DEJ2000",
            "frame": "icrs",
            "id_column": "WISE",
        },
        "sdss": {
            "mags": ["uPmag", "gPmag", "rPmag", "iPmag", "zPmag"],
            "errors": ["e_uPmag", "e_gPmag", "e_rPmag", "e_iPmag", "e_zPmag"],
            "lon_column": "RA_ICRS",
            "lat_column": "DE_ICRS",
            "frame": "icrs",
            "id_column": "objID",
        },
        "2mass": {
            "mags": ["Jmag", "Hmag", "Kmag"],
            "errors": ["e_Jmag", "e_Hmag", "e_Kmag"],
            "lon_column": "RAJ2000",
            "lat_column": "DEJ2000",
            "frame": "icrs",
            "id_column": "2MASS",
        },
        "skymapper": {
            "mags": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
            "errors": ["e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF", "e_uPSF", "e_vPSF"],
            "lon_column": "RAICRS",
            "lat_column": "DEICRS",
            "frame": "icrs",
            "id_column": "ObjectId",
        },
        "panstarrs": {
            "mags": ["gmag", "rmag", "imag", "zmag", "ymag"],
            "errors": ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"],
            "lon_column": "RAJ2000",
            "lat_column": "DEJ2000",
            "frame": "icrs",
            "id_column": "objID",
        },
    },
    "positional": {
        "rosat": {"lon_column": "RAJ2000", "lat_column": "DEJ2000", "frame": "icrs", "id_column": "Name"},
        "erosita": {"lon_column": "RA_ICRS", "lat_column": "DE_ICRS", "frame": "icrs", "id_column": "IAUName"},
    },
}
