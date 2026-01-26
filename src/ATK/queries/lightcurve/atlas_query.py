import re
import time
from io import StringIO

import pandas as pd

from ...structures.Target import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request
from .lightcurve_core import get_lightcurves


def query(target: Target, **kwargs: dict):
    """
    Performs an ATLAS light curve query

    kwargs:
    - limit_mjd: if True, limit MJD to mid 2025 for testing
    """

    url = "https://fallingstar-data.com/forcedphot"

    if not kwargs.get("username"):
        raise ValueError(
            f"An ATLAS forcedphot ({url}) username is required for ATLAS light curve queries. One should be provided with 'username=...'."
        )
    if not kwargs.get("password"):
        raise ValueError(
            f"An ATLAS forcedphot ({url}) password is required for ATLAS light curve queries. One should be provided with 'password=...'."
        )

    response = send_request(
        "atlas",
        f"{url}/api-token-auth/",
        method="POST",
        data={"username": kwargs["username"], "password": kwargs["password"]},
        message="ATLAS API authorization failed. Please check that you have provided valid login details.",
    )
    if response is RETURNS.EXCEPTION:
        return response

    token = response.json()["token"]
    headers = {"Authorization": f"Token {token}", "Accept": "application/json"}
    epoch = target.coords.obstime.datetime.year

    request_data = {
        "ra": target.coords.ra.value,
        "dec": target.coords.dec.value,
        "mjd_min": 60858.5 if kwargs.get("limit_mjd") else 50000.0,
        "radec_epoch_year": epoch,
        "use_reduced": True,
    }

    if target.correction != "none":
        request_data["propermotion_ra"] = target.coords.pm_ra_cosdec.value
        request_data["propermotion_dec"] = target.coords.pm_dec.value

    task_url = None
    while not task_url:
        response = send_request("atlas", f"{url}/queue/", method="POST", data=request_data, headers=headers)
        if response is RETURNS.EXCEPTION:
            return response

        if response.status_code == 201:
            task_url = response.json()["url"]
        elif response.status_code == 429:
            message = response.json()["detail"]

            t_sec = re.findall(r"available in (\d+) seconds", message)
            t_min = re.findall(r"available in (\d+) minutes", message)
            t_wait = int(t_sec[0]) if t_sec else int(t_min[0]) if t_min else 10
            time.sleep(t_wait)
        else:
            return RETURNS.EXCEPTION

    result_url = None
    while not result_url:
        response = send_request("atlas", task_url, headers=headers)
        if response is RETURNS.EXCEPTION:
            return response
        if response.status_code == 200:
            if response.json()["finishtimestamp"]:
                result_url = response.json()["result_url"]
            else:
                time.sleep(10)
        else:
            return RETURNS.EXCEPTION

    data = send_request("atlas", result_url, headers=headers)
    if data is RETURNS.EXCEPTION:
        return data

    df = pd.read_csv(StringIO(data.text.replace("###", "")), sep="\\s+")
    if df.empty:
        return RETURNS.NULL

    # basic filtering (old)
    """
    df = df[np.abs(df["uJy"]) > 3]
    df = df[np.abs(df["duJy"]) < 4000]
    df = df[df["dm"] > 0]
    """

    # basic filtering (as recommended by ATLAS team)
    mask = (
        (df["duJy"] < 10000)
        & (df["err"] == 0)
        & (df["x"] > 100)
        & (df["x"] < 10460)
        & (df["y"] > 100)
        & (df["y"] < 10460)
        & (df["maj"] < 5)
        & (df["maj"] > 1.6)
        & (df["min"] < 5)
        & (df["min"] > 1.6)
        & (df["apfit"] > -1)
        & (df["apfit"] < -0.1)
        & (df["mag5sig"] > 17)
        & (df["Sky"] > 17)
    )

    df = df[mask]

    # signal to noise + negative flux filtering (magnitudes derived from negative fluxes not physical)
    snr = df["uJy"] / df["duJy"]
    df = df[(df["uJy"] > 0) & (snr >= 3)]

    df = df.rename(columns={"MJD": "mjd", "m": "mag", "dm": "mag_err", "F": "band", "RA": "ra", "Dec": "dec"})

    lcs = get_lightcurves("atlas", target, kwargs["radius"], df, kwargs.get("split", False))

    return lcs
