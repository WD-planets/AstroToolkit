from enum import Enum, auto


class RETURNS(Enum):
    SUCCESS = auto()
    NULL = auto()
    EXCEPTION = auto()


GAIA_CATALOGUE = "I/355/gaiadr3"

QUERY_KINDS = ["vizier", "image", "lightcurve", "spectrum", "sed", "hrd"]

PLOT_DIMENSIONS = {"image": [2, 2]}
