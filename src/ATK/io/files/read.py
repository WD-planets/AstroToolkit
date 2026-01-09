import warnings
from pathlib import Path
from types import NoneType

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU
from astropy.table import Table
from astropy.units import Quantity, Unit
from astropy.wcs import WCS

from ...configuration.base_config import translator
from ...structures.definitions import Image, QueryResult
from ...utilities.mapping import build_structure_map, get_query_result_map
from ...utilities.misc import get_package_version

SKYCOORD_KEYS = ("ATK_RA", "ATK_DEC", "ATK_PMRA", "ATK_PMDEC", "ATK_DISTANCE", "ATK_FRAME", "ATK_EPOCH")

# ---------
# UTILITIES
# ---------


def read_quantity(header: Header, key: str, default_unit: Unit | None = None) -> Quantity:
    """
    Attempts to read an astropy quantity from a fits header (from two keys, ATK_... with unit ATK_..._U). A default unit may be provided in the case of a missing unit key
    """

    val_key = key.upper()
    unit_key = f"{key}_U"

    if val_key not in header:
        raise KeyError(f"Missing FITS keyword: {val_key}")

    if unit_key not in header:
        if not default_unit:
            raise KeyError(f"Missing FITS keyword: {unit_key}")

        if isinstance(header[val_key], (int, float)):
            return header[val_key] * default_unit
        else:
            return header[val_key]
    else:
        try:
            unit = u.Unit(header[unit_key])
        except Exception:
            raise ValueError(f"Invalid unit in key '{unit_key}': {header[unit_key]}")

    return header[val_key] * unit


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
        ra=read_quantity(hdr, "ATK_RA", u.deg),
        dec=read_quantity(hdr, "ATK_DEC", u.deg),
        pm_ra_cosdec=read_quantity(hdr, "ATK_PMRA", u.mas / u.yr),
        pm_dec=read_quantity(hdr, "ATK_PMDEC", u.mas / u.yr),
        distance=read_quantity(hdr, "ATK_DISTANCE", u.pc),
        frame=hdr.get("ATK_FRAME"),
        obstime=hdr.get("ATK_EPOCH"),
    )

    return coord


def parse_primary_header(path: str | Path, hdr: Header) -> QueryResult:
    """
    Parse the primary HDU header, which contains all query information
    """

    # warn user if there is an ATK version mismatch in case
    file_ver = hdr.get("ATK_VER")
    if not file_ver:
        warnings.warn(f"Could not find ATK_VER key in file '{path}', and so the file may not be read correctly.")
    if file_ver and file_ver != get_package_version():
        warnings.warn(f"File '{path}' was generated with a different ATK version, and so may not be read correctly.")

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
        if f"{key}_U" in ATK_keys:
            val = read_quantity(hdr, key)
        else:
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

    # get container constructor
    ctnr_constr = structure_map.get(hdr.get("ATK_KIND"))

    # if no container exists (i.e. in vizier queries), just set .data = dataframe
    if not ctnr_constr:
        return df

    # populate dict with dataframe columns as np arrays
    ctnr_data = {}
    for col in df:
        if hasattr(ctnr_constr, col):
            ctnr_data[col] = np.asarray(df[col])

    # construct container with data
    ctnr = ctnr_constr(**ctnr_data)

    # parse header keys into container's non-data attributes
    ctnr = parse_header(path, hdr, ctnr)

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

        # all extensions must contain this key
        if not hdu.header.get("ATK_EXT", None):
            warnings.warn(f"ATK: Non-ATK extension found in target file {path} has been ignored.")
            continue

        # images use ImageHDU + BinTableHDU (for overlay), everything else uses BinTableHDU
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
            overlay = parse_generic_bintable(structure, path, hdul_no_primary[index + 1].header, hdul_no_primary[index + 1].data)
            data.overlay = None if overlay.empty else overlay
            completed.append(hdul_no_primary[index + 1])

        structure.data.append(data)

    return structure
