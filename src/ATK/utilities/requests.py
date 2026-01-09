import warnings

import requests
from requests import Session
from requests.adapters import HTTPAdapter, Retry
from requests.models import Response

from ..utilities.defaults import CONNECTION_ERRORS, RETURNS


def print_bad_response(survey: str, response: Response):
    """
    Attempts to print the problem(s) encountered with a 'bad' request
    """

    print(f"Note: experiencing issues with {survey} (bad status code {response.status_code}).")

    try:
        for key, val in response.json().items():
            print(f"    {key}:", val)
        print()
        return
    except ValueError:
        pass

    try:
        print(f"    {response.text}\n")
        return
    except Exception:
        return


def send_request(survey: str, url: str, method: str = "GET", message=None, **kwargs) -> Response | RETURNS:
    """
    Fetches data from a given URL. Supports HTTP verbs via method
    """

    s = Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods={"GET", "POST"})
    # s.mount("http://", HTTPAdapter(max_retries=retries))
    s.mount("https://", HTTPAdapter(max_retries=retries))

    # send request
    try:
        response = s.request(method.upper(), url, timeout=180, **kwargs)
        response.raise_for_status()

    except CONNECTION_ERRORS as e:
        print(e)
        if message:
            print(message)
        else:
            warnings.warn(f"Note: experiencing issues with {survey} (timeout)")
        return RETURNS.EXCEPTION

    except requests.HTTPError:
        if message:
            print(message)
        else:
            print_bad_response(survey, response)

        return RETURNS.EXCEPTION

    return response
