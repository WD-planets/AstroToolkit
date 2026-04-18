from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_overlays.ini"

DEFAULTS = {
    "photometric": {
        "gaia": {
            "mags": ["Gmag", "BPmag", "RPmag"],
            "errors": ["e_Gmag", "e_BPmag", "e_RPmag"],
            "lon": "RA_ICRS",
            "lat": "DE_ICRS",
            "frame": "icrs",
        },
        "galex": {
            "mags": ["NUVmag", "FUVmag"],
            "errors": ["e_NUVmag", "e_FUVmag"],
            "lon": "RAJ2000",
            "lat": "DEJ2000",
            "frame": "icrs",
        },
        "wise": {
            "mags": ["W1mag", "W2mag", "W3mag", "W4mag"],
            "errors": ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"],
            "lon": "RAJ2000",
            "lat": "DEJ2000",
            "frame": "icrs",
        },
        "sdss": {
            "mags": ["uPmag", "gPmag", "rPmag", "iPmag", "zPmag"],
            "errors": ["e_uPmag", "e_gPmag", "e_rPmag", "e_iPmag", "e_zPmag"],
            "lon": "RA_ICRS",
            "lat": "DE_ICRS",
            "frame": "icrs",
        },
        "2mass": {
            "mags": ["Jmag", "Hmag", "Kmag"],
            "errors": ["e_Jmag", "e_Hmag", "e_Kmag"],
            "lon": "RAJ2000",
            "lat": "DEJ2000",
            "frame": "icrs",
        },
        "skymapper": {
            "mags": ["gPSF", "rPSF", "iPSF", "zPSF", "uPSF", "vPSF"],
            "errors": ["e_gPSF", "e_rPSF", "e_iPSF", "e_zPSF", "e_uPSF", "e_vPSF"],
            "lon": "RAICRS",
            "lat": "DEICRS",
            "frame": "icrs",
        },
        "panstarrs": {
            "mags": ["gmag", "rmag", "imag", "zmag", "ymag"],
            "errors": ["e_gmag", "e_rmag", "e_imag", "e_zmag", "e_ymag"],
            "lon": "RAJ2000",
            "lat": "DEJ2000",
            "frame": "icrs",
        },
    },
    "positional": {
        "rosat": {"lon": "RAJ2000", "lat": "DEJ2000", "frame": "icrs"},
        "erosita": {"lon": "RA_ICRS", "lat": "DE_ICRS", "frame": "icrs"},
    },
}
