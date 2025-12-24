import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU
from astropy.table import Table
from astropy.units import Quantity, UnitBase
from astropy.wcs import WCS

from ...configuration.base_config import translator
from ...structures.definitions import Image, QueryResult
from ...utilities.mapping import build_structure_map, get_query_result_map

SKYCOORD_KEYS = ("ATK_RA", "ATK_DEC", "ATK_PMRA", "ATK_PMDEC", "ATK_DISTANCE", "ATK_FRAME", "ATK_EPOCH")

# ---------
# UTILITIES
# ---------


def get_header_quantity(hdr: Header, key: str, unit: UnitBase) -> Quantity | None:
    """
    Safely fetch a value which has an intended astropy unit from a fits header
    """

    val = hdr.get(key)
    if not val:
        return val

    return val * unit


# -----------------
# PARSING FUNCTIONS
# -----------------


def parse_header_skycoord(path: str | Path, hdr: Header) -> SkyCoord:
    """
    Parses multiple header keys into an astropy SkyCoord object
    """

    # check for missing SkyCoord keys
    if any([key not in hdr.keys() for key in SKYCOORD_KEYS]):
        raise ValueError(f"Invalid SkyCoord decomposition in header of {path}")

    coord = SkyCoord(
        ra=get_header_quantity(hdr, "ATK_RA", u.deg),
        dec=get_header_quantity(hdr, "ATK_DEC", u.deg),
        pm_ra_cosdec=get_header_quantity(hdr, "ATK_PMRA", u.mas / u.yr),
        pm_dec=get_header_quantity(hdr, "ATK_PMDEC", u.mas / u.yr),
        distance=get_header_quantity(hdr, "ATK_DISTANCE", u.pc),
        frame=hdr.get("ATK_FRAME"),
        obstime=hdr.get("ATK_EPOCH"),
    )

    return coord


def parse_primary_header(path: str | Path, hdr: Header) -> QueryResult:
    """
    Parse the primary HDU header, which contains all query information
    """

    structure = get_query_result_map()[hdr.get("ATK_KIND")]()

    # a position should always be generated, so this should not trigger
    if "ATK_SC" not in hdr:
        raise ValueError(f"No targeting information found in primary header of {path}.")

    structure.position = parse_header_skycoord(path, hdr)

    structure = parse_header(path, hdr, structure)

    return structure


def parse_header(path: str | Path, hdr: Header, obj: object) -> object:
    """
    Parse fits header keys into their corresponding object attributes
    """

    sc_attr = hdr.get("ATK_SC", None)
    if sc_attr:
        setattr(obj, sc_attr, parse_header_skycoord(path, hdr))

    ATK_keys = {key: val for key, val in hdr.items() if key.startswith("ATK_")}
    for key, val in ATK_keys.items():
        # use same translator as base config to get correct data types from ATK header keys
        val = translator(val)

        # get rid of ATK_ prefix and lower key to match attribute
        attr = key[4:].lower()

        if hasattr(obj, attr):
            setattr(obj, attr, val)

    return obj


def parse_generic_bintable(structure: QueryResult, path: str | Path, hdr: Header, data: any) -> any:
    """
    Parse a bintable extension into either a data container or dataframe (for vizier queries)
    """

    structure_map = build_structure_map()
    df = Table(data).to_pandas()
    ctnr = structure_map.get(hdr.get("ATK_KIND"), lambda: None)()

    # if no container exists (i.e. in vizier queries), just set .data = dataframe
    if not ctnr:
        return df

    # otherwise, parse header keys into container
    ctnr = parse_header(path, hdr, ctnr)

    # populate container with dataframe columns as np arrays
    for col in df:
        if hasattr(ctnr, col):
            setattr(ctnr, col, np.asarray(df[col]))

    return ctnr


def parse_generic_imagehdu(structure: QueryResult, path: str | Path, hdu: ImageHDU) -> any:
    """
    Parse an imageHDU into a data container (i.e. an Image)
    """

    structure_map = build_structure_map()
    ctnr = structure_map.get(hdu.header.get("ATK_KIND"), lambda: None)()

    if not ctnr or not isinstance(ctnr, Image):
        raise Exception(f"Could not parse ImageHDU in {path} into Image.")

    ctnr = parse_header(path, hdu.header, ctnr)

    ctnr.hdu = hdu
    ctnr.wcs = WCS(hdu.header)

    return ctnr


# ----
# MAIN
# ----


def read_local(path: str | Path) -> QueryResult:
    """
    Read a local ATK fits file back into the original data structure that was used to generate it
    """

    hdul = fits.open(path)

    structure = parse_primary_header(path, hdul[0].header)

    # iterate through extensions
    completed = []
    hdul_no_primary = hdul[1:]
    for index, hdu in enumerate(hdul_no_primary):
        if hdu in completed:
            continue

        if not hdu.header.get("ATK_EXT", None):
            warnings.warn(f"ATK: Non-ATK extension found in target file {path} has been ignored.")
            continue
        if not isinstance(hdu, (BinTableHDU, ImageHDU)):
            warnings.warn(f"ATK: Unexpected table type '{type(hdu)}' in target file {path} has been ignored.")
            continue

        # everything except images
        if isinstance(hdu, BinTableHDU):
            data = parse_generic_bintable(structure, path, hdu.header, hdu.data)
            completed.append(hdu)

        # images
        elif isinstance(hdu, ImageHDU):
            # get image hdu
            data = parse_generic_imagehdu(structure, path, hdu)
            completed.append(hdu)
            # get overlay from next extension
            data.overlay = parse_generic_bintable(structure, path, hdul_no_primary[index + 1].header, hdul_no_primary[index + 1].data)
            completed.append(hdul_no_primary[index + 1])

        structure.data.append(data)

    return structure
