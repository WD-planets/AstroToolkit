import typing
from dataclasses import fields
from types import NoneType, UnionType
from typing import get_args, get_origin

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU, PrimaryHDU
from astropy.table import Table
from astropy.wcs import WCS

from .definitions import Image

# types (in typehints) that should be considered as being columns of a dataframe
COLUMN_TYPES = (np.ndarray, pd.Series, list, tuple, set)

BASIC_TYPES = (int, float, str, bool, NoneType)

# ----------------------
# HEADER WRITE FUNCTIONS
# ----------------------


def write_fallback(attr: str, hdr: Header, key: str, value: any) -> Header:
    """
    Writes generic data types to the header (e.g. str, int, float)
    """

    try:
        hdr.append((f"ATK_{key.upper()}", value))
    except Exception:
        pass

    try:
        hdr.append((f"ATK_{key.upper()}", str(value)))
    except Exception:
        raise ValueError(f"Failed to write value '{value}' of type '{type(value)}' in attribute '{attr}' to FITS header key '{key}'.")

    return hdr


def write_skycoord(attr: str, hdr: Header, coord: SkyCoord) -> Header:
    """
    Splits a SkyCoord into its components and writes this as a set of header keys
    """

    hdr.append(("ATK_SC", attr, "Origin attribute of decomposed SkyCoord"))
    hdr.append(("ATK_RA", coord.ra.value, "Source RA (deg)"))
    hdr.append(("ATK_DEC", coord.dec.value, "Source DEC (deg)"))
    hdr.append(("ATK_FRAME", coord.frame.name, "Coordinate frame of ATK_RA and ATK_DEC"))

    if coord.obstime:
        hdr.append(("ATK_EPOCH", coord.obstime.fits, "Epoch of ATK_RA and ATK_DEC"))

    # proper motion data
    if coord.data.differentials:
        hdr.append(("ATK_PMRA", coord.pm_ra_cosdec.value, "Source Proper Motion in RA (mas/yr)"))
        hdr.append(("ATK_PMDEC", coord.pm_dec.value, "Source Proper Motion in DEC (mas/yr)"))
    else:
        hdr.append(("ATK_PMRA", None))
        hdr.append(("ATK_PMDEC", None))

    # distance
    if coord.distance != u.one:
        hdr.append(("ATK_DISTANCE", coord.distance.value, "Source distance (1/p) (pc)"))
    else:
        hdr.append(("ATK_DISTANCE", None))

    return hdr


# Map to special header writing functions
WRITE_MAP = {SkyCoord: lambda **kwargs: write_skycoord(kwargs["attr"], kwargs["hdr"], kwargs["value"])}

# ---------
# UTILITIES
# ---------


def is_strict_column_type(typ: UnionType | type) -> bool:
    """
    Returns True if all non-None types in a typehint are array-like
    """

    # Unions (i.e. x | y | z)
    if isinstance(typ, UnionType):
        args = [a for a in get_args(typ) if a is not type(None)]
        if not args:
            return False
        return all(is_strict_column_type(a) for a in args)

    # Generic nested types (i.e. list[float])
    origin = get_origin(typ)
    if origin in COLUMN_TYPES:
        return True

    # Basic non-nested type (i.e. np.ndarray)
    return typ in COLUMN_TYPES


def get_cols(structure: any) -> tuple[str]:
    """
    Returns a list of array-like attributes of a data structure using typehinting
    """

    cols = []

    hints = typing.get_type_hints(structure.__class__)

    for field in fields(structure):
        field_type = hints.get(field.name)
        if field_type is None:
            continue

        if is_strict_column_type(field_type):
            cols.append(field.name)

    return cols


# -----------------------
# GENERAL TRANSFORMATIONS
# -----------------------


def struct_to_dataframe(structure: any) -> pd.DataFrame:
    """
    Combines the array-like attributes of a data structure into a single pandas DataFrame
    """

    cols = get_cols(structure)

    data = {}
    for col in cols:
        val = getattr(structure, col)
        if val is None:
            data[col] = np.empty(0)
            continue

        if not isinstance(val, COLUMN_TYPES):
            val = [val]
        data[col] = val

    return pd.DataFrame.from_dict(data)


def struct_to_hdu(structure: any, ignore_attrs: list = [], hdu_kind: PrimaryHDU | BinTableHDU | ImageHDU = BinTableHDU) -> BinTableHDU:
    """
    Convert a data structure into a fits HDU.
    """

    # make sure kind is always ignored by dispatcher, this is handled separately below
    ignore_attrs.append("kind")

    hdr = Header()

    # tag all HDUs as coming from ATK
    hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))

    # combine array-like attributes into a dataframe
    cols = get_cols(structure)
    df = struct_to_dataframe(structure)
    tbl = Table.from_pandas(df)

    # PrimaryHDU stores query kind, extensions store data container kind
    kind_str = "ATK query kind" if hdu_kind is PrimaryHDU else "ATK container kind"
    kind = structure.__dict__.get("kind", type(structure).__name__)
    hdr.append(("ATK_KIND", kind, kind_str))

    # iterate through remaining structure attributes and write them to the header
    for attr, val in structure.__dict__.items():
        if attr in ignore_attrs:
            continue
        if attr in cols:
            continue

        hdr = WRITE_MAP.get(type(val), write_fallback)(attr=attr, hdr=hdr, key=attr, value=val)

    # generate HDU
    if hdu_kind == PrimaryHDU:
        hdu = PrimaryHDU(None, header=hdr)
    elif hdu_kind == ImageHDU:
        hdu = ImageHDU(tbl, header=hdr, name=structure.__str__())
    elif hdu_kind == BinTableHDU:
        hdu = BinTableHDU(tbl, header=hdr, name=structure.__str__())
    else:
        raise ValueError(f"Invalid table type '{kind}' passed to struct_to_hdu.")

    return hdu


# -----------------------
# SPECIAL TRANSFORMATIONS
# -----------------------


def image_to_hdu(image: Image):
    hdu = image.hdu
    hdr = hdu.header

    hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))
    hdr.append(("ATK_KIND", "Image", "ATK container kind"))
    for attr, val in image.__dict__.items():
        # basic types and those with special writing functions (e.g. SkyCoord decomposition)
        if type(val) in WRITE_MAP or type(val) in BASIC_TYPES:
            hdr = WRITE_MAP.get(type(val), write_fallback)(attr=attr, hdr=hdr, key=attr, value=val)

    image_hdu = ImageHDU(data=hdu.data, header=hdr, name=image.__str__())

    overlay_hdr = Header()
    overlay_hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))
    overlay = image.overlay
    if overlay is None:
        overlay = pd.DataFrame()
    table = Table.from_pandas(overlay)
    overlay_hdu = BinTableHDU(data=table, header=overlay_hdr, name="<Overlay Data>")

    return (image_hdu, overlay_hdu)
