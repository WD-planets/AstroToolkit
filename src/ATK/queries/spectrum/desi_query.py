import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astroquery.utils.tap.core import TapPlus
from sparcl.client import SparclClient

from ...structures.definitions import Spectrum, Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from ...utilities.misc import angle_to_quantity, suppress_stdout


def desi_cone_search(position: SkyCoord, radius: Quantity, table="desi_dr1.zpix", columns="glon, glat, targetid"):
    """
    Performs a Tap query to DESI, returning any objects with spectra within a cone
    """

    tap = TapPlus(url="https://datalab.noirlab.edu/tap")
    tap.TIMEOUT = 180

    radius_deg = radius.to(u.deg).value

    # DESI table only has galactic coords
    position = position.galactic

    # query string, uses trig to make a cone search since this is not supported by DESI by default
    # only returns primary spectra with no warning flag
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

    # launch query
    try:
        job = tap.launch_job(query)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    # get results and convert to DataFrame
    result = job.get_results()
    df = result.to_pandas()

    if df.empty:
        return RETURNS.NULL

    # compute separations explicitly (no clipping)
    result_coords = SkyCoord(l=df["glon"].to_numpy() * u.deg, b=df["glat"].to_numpy() * u.deg, frame="galactic")

    # get separation from cone centre in arcsec + sort
    df["sep_arcsec"] = position.separation(result_coords).to(u.arcsec).value
    df = df.sort_values("sep_arcsec").reset_index(drop=True)

    return df


def retrieve_desi_spectra(targetids: list):
    """
    Retrieves spectra for a list of DESI target IDs
    """

    # suppresses any announcements not covered by verbose=False
    with suppress_stdout():
        client = SparclClient(verbose=False)

    # get DESI fields
    include = client.get_all_fields(dataset_list=["DESI-DR1"])

    # get spectra
    response = client.retrieve_by_specid(specid_list=targetids, dataset_list=["DESI-DR1"], include=include)

    return response


def query(target: Target, **kwargs: dict):
    """
    Performs a DESI spectrum query
    """

    # get detections by cone search
    data = desi_cone_search(target.coords, kwargs["radius"])

    # if no data is returned or an exception is encountered
    if not isinstance(data, pd.DataFrame):
        return data

    # use SparclClient to retreive spectra for detections
    records = retrieve_desi_spectra(data["targetid"].tolist())
    if not records[0]["status"]["success"]:
        return RETURNS.EXCEPTION

    # get necessary information from spectra and construct array of Spectrum objects
    spectra = []
    for record in records.data[1:]:
        if record["specprimary"]:
            spec = Spectrum(
                "desi",
                wavelength=record["wavelength"] * u.Unit("Angstrom"),
                flux=record["flux"] * u.Unit("1e-17 erg cm-2 s-1 Angstrom-1"),
                exposure=record["exptime"] * u.s,
            )
            spec_pos = SkyCoord(ra=record["ra"] * u.deg, dec=record["dec"] * u.deg, frame="icrs")
            spec.position = spec_pos
            spec.separation = angle_to_quantity(spec_pos.separation(target.coords), kwargs["radius"].unit)
            spec.program = record["program"]
            spectra.append(spec)

    return spectra
