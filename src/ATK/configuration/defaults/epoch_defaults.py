from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_epochs.ini"

DEFAULTS = {
    "vizier_aliases": {
        "gaia": "2016-01-01T00:00:00.000",
        "panstarrs": "2012-01-01T00:00:00.000",
        "skymapper": "2016-01-01T00:00:00.000",
        "galex": "2006-08-01T00:00:00.000",
        "rosat": "1991-01-01T00:00:00.000",
        "sdss": "2006-01-01T00:00:00.000",
        "wise": "2010-06-01T00:00:00.000",
        "2mass": "1999-01-01T00:00:00.000",
        "erosita": "2022-01-01T00:00:00.000",
    },
    "lightcurve_surveys": {
        "ztf": "2019-01-01T00:00:00.000",
        "atlas": "2021-01-01T00:00:00.000",
        "gaia_lc": "2016-01-01T00:00:00.000",
        "asassn": "2015-01-01T00:00:00.000",
        "crts": "2008-01-01T00:00:00.000",
        "tess": "2020-01-01T00:00:00.000",
    },
    "spectrum_surveys": {"sdss": "2017-01-01T00:00:00.000"},
}
