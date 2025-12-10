from requests import Session
from requests.adapters import HTTPAdapter, Retry
from requests.models import Response

from ..utilities.defaults import RETURNS


def send_request(survey: str, url: str, method: str = "GET", **kwargs) -> Response | None:
    """
    Fetches data from a given URL with retry logic.
    Supports GET and POST (or other HTTP verbs via `method`).
    """

    s = Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    s.mount("http://", HTTPAdapter(max_retries=retries))
    s.mount("https://", HTTPAdapter(max_retries=retries))

    try:
        response = s.request(method.upper(), url, timeout=180, **kwargs)
    except TimeoutError:
        print(f"Note: experiencing issues with {survey} (timeout)")
        return RETURNS.EXCEPTION

    if response.status_code != 200:
        print(f"Note: experiencing issues with {survey} (bad status code {response.status_code})")
        return RETURNS.EXCEPTION

    return response
