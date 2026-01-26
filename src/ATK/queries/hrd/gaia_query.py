import astropy.units as u
import numpy as np

from ...structures.HRD import HRD
from ...structures.Target import Target
from ...Tools.query import query as general_query
from ...utilities.defaults import RETURNS


def query(target: Target, **kwargs):
    bands = kwargs["colour"].split("-")
    for band in bands:
        if band not in ["Gmag", "BPmag", "RPmag"]:
            raise ValueError(f"Unknown Gaia band '{band}'.")
    if not target.identifier:
        raise ValueError("Targets must be a source ID or list of source IDs.")

    gaia_data = general_query("vizier", target=target, survey="gaia")
    if not gaia_data.data or gaia_data.exception:
        return RETURNS.EXCEPTION

    data = gaia_data.data[0].data
    plx = data["Plx"][0]
    distance = 1 / (plx * 1e-3)
    if np.isnan(plx):
        return RETURNS.NULL

    abs_mag = data[kwargs["mag"]][0] + 5 * np.log10(plx / 1000) + 5
    colour = data[bands[0]][0] - data[bands[1]][0]

    hrd = HRD(
        survey="gaia",
        abs_mag_band=kwargs["mag"],
        colour_bands=kwargs["colour"],
        colour=np.asarray([colour]),
        distance=np.asarray([distance]) * u.pc,
        abs_mag=np.asarray([abs_mag]),
        correction="n/a",
        identifier=target.identifier,
    )

    return [hrd]
