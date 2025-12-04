import typing
from dataclasses import fields

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.table import Table
from astropy.wcs import WCS

# types (in typehints) that should be considered as being columns of a dataframe
COLUMN_TYPES = (np.ndarray, pd.Series, list)

# ----------------------
# HEADER WRITE FUNCTIONS
# ----------------------


def write_fallback(hdr: Header, key: str, value: any):
    try:
        hdr.append((f"ATK_{key.upper()}", value))
    except Exception:
        raise ValueError(f"Failed to write value '{value}' of type '{type(value)}' to FITS header key '{key}'.")

    return hdr


def write_skycoord(hdr: Header, coord: SkyCoord):
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


def write_wcs(hdr: Header):
    return hdr


WRITE_MAP = {SkyCoord: lambda hdr, key, val: write_skycoord(hdr, val), WCS: lambda hdr, key, val: write_wcs(hdr)}

# ---------
# UTILITIES
# ---------


def get_cols(structure: any):
    cols = []

    hints = typing.get_type_hints(structure.__class__)
    for field in fields(structure):
        field_type = hints.get(field.name, None)
        if not field_type:
            continue

        args = typing.get_args(field_type)
        if field_type in COLUMN_TYPES or any(COL in args for COL in COLUMN_TYPES):
            cols.append(field.name)

    return cols


# -----------------------
# GENERAL TRANSFORMATIONS
# -----------------------


def struct_to_dataframe(structure: any):
    cols = get_cols(structure)

    data = {}
    for col in cols:
        val = getattr(structure, col)
        if not val:
            data[col] = np.empty(0)
            continue

        if not isinstance(val, COLUMN_TYPES):
            val = [val]
        data[col] = val

    return pd.DataFrame.from_dict(data)


def struct_to_hdu(structure: any, ignore_attrs: tuple = (), kind: BinTableHDU | PrimaryHDU = BinTableHDU) -> BinTableHDU:
    hdr = Header()

    cols = get_cols(structure)
    df = struct_to_dataframe(structure)
    tbl = Table.from_pandas(df)

    for attr, val in structure.__dict__.items():
        if attr in ignore_attrs:
            continue
        if attr in cols:
            continue

        hdr = WRITE_MAP.get(type(val), write_fallback)(hdr, attr, val)

    if kind == PrimaryHDU:
        hdu = PrimaryHDU(tbl, header=hdr)
    elif kind == BinTableHDU:
        hdu = BinTableHDU(tbl, header=hdr, name=structure.__str__())
    else:
        raise ValueError("Invalid table type '{kind}' passed to struct_to_hdu.")

    return hdu
