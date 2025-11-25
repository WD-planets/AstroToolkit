from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_aliases.ini"

DEFAULTS = {
    "vizier_aliases": {
        "gaia": "I/355/gaiadr3",
        "panstarrs": "II/349/ps1",
        "skymapper": "II/379/smssdr4",
        "galex": "II/335/galex_ais",
        "rosat": "IX/11/rosatsrc",
        "sdss": "V/154/sdss16",
        "wise": "II/311/wise",
        "2mass": "II/246/out",
        "erosita": "J/A+A/682/A34/erass1-m",
    }
}
