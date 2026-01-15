import astropy.units as u
import numpy as np

from ...structures.definitions import HRD, Target
from ...Tools import query as general_query
from ...utilities.defaults import RETURNS


def query(target: Target, **kwargs):
    bands = kwargs["colour"].split("-")
    for band in bands:
        if band not in ["g", "bp", "rp"]:
            raise ValueError(f"Unknown Gaia band '{band}'.")
    if not target.identifier:
        raise ValueError("Target must be a source ID or list of source IDs.")

    gaia_data = general_query("vizier", target=target, survey="gaia")
    if not gaia_data.data or gaia_data.exception:
        return RETURNS.EXCEPTION

    data = gaia_data.data[0]
    plx = data["Plx"][0]
    distance = 1 / (plx * 1e-3)
    if np.isnan(plx):
        return RETURNS.NULL
    abs_mag = data[f"{kwargs['mag'].upper()}mag"][0] + 5 * np.log10(plx / 1000) + 5

    colour = data[f"{bands[0].upper()}mag"][0] - data[f"{bands[1].upper()}mag"][0]

    return HRD(
        target._key,
        "gaia",
        kwargs["mag"].upper(),
        kwargs["colour"],
        np.asarray([colour]),
        np.asarray([distance]) * u.pc,
        np.asarray([abs_mag]),
    )
