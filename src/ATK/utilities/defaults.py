from enum import Enum, auto


class RETURNS(Enum):
    SUCCESS = auto()
    NULL = auto()
    EXCEPTION = auto()


GAIA_CATALOGUE = "I/355/gaiadr3"

QUERY_KINDS = ["vizier", "image", "lightcurve", "spectrum", "sed", "hrd"]

PLOT_DIMENSIONS = {"image": [2, 2]}

# this will need updating, also haven't checked
EFFECTIVE_WAVELENGTHS = {
    "panstarrs": {"g": 481.0, "r": 615.5, "i": 750.3, "z": 866.8, "y": 961.4},
    "galex": {"FUV": 154.9, "NUV": 230.3},
    "sdss": {"u": 355.1, "g": 468.6, "r": 616.6, "i": 748.0, "z": 893.2},
    "2mass": {"J": 1069.1, "H": 1446.5, "Ks": 2155.8},
    "wise": {"W1": 3400.0, "W2": 4600.0, "W3": 12000.0, "W4": 22000.0},
    "dss": {"Blue": 395.0, "Red": 630.0, "IR": 800.0},
    "skymapper": {"u": 350.0, "v": 450.0, "g": 480.0, "r": 625.0, "i": 775.0, "z": 870.0},
}
