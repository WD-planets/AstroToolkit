from enum import Enum, auto
from urllib.error import HTTPError

from requests.exceptions import ConnectionError, ConnectTimeout
from requests.exceptions import HTTPError as requests_HTTPError
from requests.exceptions import RetryError


class RETURNS(Enum):
    NULL = auto()
    EXCEPTION = auto()


GAIA_CATALOGUE = "I/355/gaiadr3"

QUERY_KINDS = ["vizier", "image", "lightcurve", "spectrum", "sed", "hrd"]

PLOT_DIMENSIONS = {"image": [2, 2], "sed": [3, 2], "spectrum": [2, 1], "lightcurve": [2, 1], "hrd": [2, 2]}

CONNECTION_ERRORS = (TimeoutError, ConnectionError, ConnectTimeout, HTTPError, requests_HTTPError, RetryError)
