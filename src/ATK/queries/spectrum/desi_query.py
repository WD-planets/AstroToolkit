import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.utils.tap.core import TapPlus
from sparcl.client import SparclClient

from ...structures.definitions import Spectrum, Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from ...utilities.misc import suppress_stdout


def desi_cone_search(position: SkyCoord, radius: float, table="desi_dr1.zpix", columns="glon, glat, targetid"):
    tap = TapPlus(url="https://datalab.noirlab.edu/tap")
    tap.TIMEOUT = 180

    radius_deg = radius / 3600.0
    position = position.galactic

    query = f"""
        SELECT {columns}
        FROM {table} AS t
        WHERE 57.29577951308232 * ACOS(
            SIN(RADIANS({position.b.value})) * SIN(RADIANS(glat)) +
            COS(RADIANS({position.b.value})) * COS(RADIANS(glat)) *
            COS(RADIANS(glon - {position.l.value}))
        ) < {radius_deg}
        AND t.survey = 'main'
        AND t.zwarn = 0
        AND t.zcat_primary
    """

    query = f"""
        SELECT {columns}
        FROM {table} AS t
        WHERE 57.29577951308232 * ACOS(
            SIN(RADIANS({position.b.value})) * SIN(RADIANS(glat)) +
            COS(RADIANS({position.b.value})) * COS(RADIANS(glat)) *
            COS(RADIANS(glon - {position.l.value}))
        ) < {radius_deg}
        AND t.zwarn = 0
        AND t.zcat_primary = 'True'
    """

    # ValueError can occur if no data is returned
    try:
        job = tap.launch_job(query)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    result = job.get_results()
    df = result.to_pandas()

    # compute separations explicitly (no clipping)
    result_coords = SkyCoord(l=df["glon"].to_numpy() * u.deg, b=df["glat"].to_numpy() * u.deg, frame="galactic")

    df["sep_arcsec"] = position.separation(result_coords).to(u.arcsec).value
    df = df.sort_values("sep_arcsec").reset_index(drop=True)

    return df


def retrieve_desi_spectra(targetids: list):
    with suppress_stdout():
        client = SparclClient(verbose=False)
    include = client.get_all_fields(dataset_list=["DESI-DR1"])
    response = client.retrieve_by_specid(specid_list=targetids, dataset_list=["DESI-DR1"], include=include)

    return response


def query(target: Target, **kwargs: dict):
    data = desi_cone_search(target.coords, kwargs["radius"])

    if not isinstance(data, pd.DataFrame):
        return data
    if data.empty:
        return RETURNS.NULL

    records = retrieve_desi_spectra(data["targetid"].tolist())
    if not records[0]["status"]["success"]:
        return RETURNS.EXCEPTION

    spectra = []
    for record in records.data[1:]:
        if record["specprimary"]:
            spec = Spectrum("desi", record["wavelength"], record["flux"], record["exptime"])
            spec.program = record["program"]
            spectra.append(spec)

    return spectra
