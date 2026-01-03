from io import BytesIO

import pandas as pd

from ...structures.definitions import Target
from ...utilities.defaults import RETURNS
from ...utilities.requests import send_request


def query(target: Target, **kwargs: dict):
    radius = kwargs.get("radius") / 3600
    url = f"https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE {target.coords.ra.value} {target.coords.dec.value} {radius}&BANDNAME=g,r,i&FORMAT=CSV"

    response = send_request("ztf", url)
    if response is RETURNS.EXCEPTION:
        return response

    data = pd.read_csv(BytesIO(response.content))
    if not len(data):
        return RETURNS.NULL

    for oid in set(data["oid"].tolist()):
        obj_data = data[data["oid"] == oid]
