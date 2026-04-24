import time
from io import StringIO
from pathlib import Path

import astropy.units as u
import pandas as pd
from bs4 import BeautifulSoup

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request
from .lightcurve_core import get_lightcurves

# query rate limit in seconds
CRTS_MIN_INTERVAL = 15
TIMER_PATH = Path.home() / ".AstroToolkit" / "ATK_CRTS_TIMER.ini"

CRTS_BASE_URL = "http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv_new.cgi"


def rate_limit():
    """
    Enforce CRTS 15s query limit using a timer file
    """

    now = time.time()

    if TIMER_PATH.exists():
        last = float(TIMER_PATH.read_text())
        elapsed = now - last
        if elapsed < CRTS_MIN_INTERVAL:
            wait = CRTS_MIN_INTERVAL - elapsed
            time.sleep(wait)

    TIMER_PATH.write_text(str(time.time()))


def fetch_csv_link(url: str) -> str | RETURNS:
    """
    Fetch CRTS HTML page and extract CSV download link
    """

    response = send_request("crts", url)
    if response is RETURNS.EXCEPTION:
        return response

    soup = BeautifulSoup(response.content, "html.parser")

    # safely grab csv link
    link = soup.find("a", href=lambda h: h and h.endswith(".csv"))

    return link["href"] if link else RETURNS.NULL


def query(target: Target, **kwargs: dict):
    """
    Perform a CRTS light curve query
    """

    radius = kwargs["radius"].to(u.arcmin).value

    rate_limit()

    url = f"{CRTS_BASE_URL}?RADec={target.coords.ra.value} {target.coords.dec.value}&Rad={radius}&OUT=csv&SHORT=short&DB=photcat"

    csv_link = fetch_csv_link(url)
    if csv_link in (RETURNS.NULL, RETURNS.EXCEPTION):
        return csv_link

    csv_response = send_request("crts", csv_link)
    if csv_response is RETURNS.EXCEPTION:
        return csv_response

    df = pd.read_csv(StringIO(csv_response.text))
    if df.empty:
        return RETURNS.NULL

    df = df.rename(columns={"MasterID": "id", "Mag": "mag", "Magerr": "mag_err", "RA": "ra", "Dec": "dec", "MJD": "mjd"})
    df["band"] = ["v"] * len(df)

    lcs = get_lightcurves("crts", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
