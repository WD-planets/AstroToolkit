import typing
from dataclasses import fields
from types import UnionType
from typing import get_args, get_origin

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU, PrimaryHDU
from astropy.table import Table
from astropy.wcs import WCS

# types (in typehints) that should be considered as being columns of a dataframe
COLUMN_TYPES = (np.ndarray, pd.Series, list, tuple, set)

# ----------------------
# HEADER WRITE FUNCTIONS
# ----------------------


def write_fallback(hdr: Header, key: str, value: any) -> Header:
    try:
        hdr.append((f"ATK_{key.upper()}", value))
    except Exception:
        raise ValueError(f"Failed to write value '{value}' of type '{type(value)}' to FITS header key '{key}'.")

    return hdr


def write_skycoord(hdr: Header, coord: SkyCoord) -> Header:
    hdr.append(("ATK_SKYCOORD", True, "If True, a SkyCoord is present in header"))
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

    return hdr


WRITE_MAP = {
    SkyCoord: lambda **kwargs: write_skycoord(kwargs["hdr"], kwargs["value"]),
    WCS: lambda **kwargs: kwargs["hdr"],
}

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


def struct_to_hdu(
    structure: any, ignore_attrs: list = [], hdu_kind: BinTableHDU | ImageHDU = BinTableHDU
) -> BinTableHDU:
    ignore_attrs.append("kind")

    hdr = Header()
    hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))

    cols = get_cols(structure)
    df = struct_to_dataframe(structure)
    tbl = Table.from_pandas(df)

    kind_str = "ATK query kind" if hdu_kind is ImageHDU else "ATK container kind"
    kind = structure.__dict__.get("kind", type(structure).__name__)
    hdr.append(("ATK_KIND", kind, kind_str))

    for attr, val in structure.__dict__.items():
        if attr in ignore_attrs:
            continue
        if attr in cols:
            continue

        hdr = WRITE_MAP.get(type(val), write_fallback)(hdr=hdr, key=attr, value=val)

    if hdu_kind == PrimaryHDU:
        hdu = PrimaryHDU(None, header=hdr)
    elif hdu_kind == ImageHDU:
        hdu = ImageHDU(tbl, header=hdr, name=structure.__str__())
    elif hdu_kind == BinTableHDU:
        hdu = BinTableHDU(tbl, header=hdr, name=structure.__str__())
    else:
        raise ValueError(f"Invalid table type '{kind}' passed to struct_to_hdu.")

    return hdu
