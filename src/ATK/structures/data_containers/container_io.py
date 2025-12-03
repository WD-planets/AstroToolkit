import typing
from dataclasses import fields

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU
from astropy.table import Table

# types (in typehints) that should be considered as being columns of a dataframe
COLUMN_TYPES = (np.ndarray, pd.Series, list)

# ----------------------
# HEADER WRITE FUNCTIONS
# ----------------------


def write_fallback(hdr: Header, key: str, value: any):
    try:
        hdr[key] = value
    except Exception:
        raise ValueError(f"Failed to write value '{value}' of type '{type(value)}' to FITS header key '{key}'.")

    return hdr


def write_skycoord(hdr: Header, coord: SkyCoord):
    print(hdr.__dict__)

    hdr["ATK_RA"] = coord.ra
    hdr["ATK_DEC"] = coord.dec

    # proper motion data
    if coord.data.differentials:
        hdr["ATK_PMRA"] = coord.pm_ra_cosdec
        hdr["ATK_PMDEC"] = coord.pm_dec
    else:
        hdr["ATK_PMRA"] = None
        hdr["ATK_PMDEC"] = None

    # distance
    if coord.distance != u.one:
        hdr["ATK_DISTANCE"] = coord.distance

    return hdr


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


def container_to_dataframe(structure: any):
    cols = get_cols(structure)

    data = {}
    for col in cols:
        val = getattr(structure, col)
        if not isinstance(val, COLUMN_TYPES):
            val = [val]
        data[col] = val

    return pd.DataFrame.from_dict(data)


def container_to_hdu(structure: any) -> BinTableHDU:
    hdr = Header()

    cols = get_cols(structure)
    df = container_to_dataframe(structure)
    tbl = Table.from_pandas(df)

    for attr, val in structure.__dict__.items():
        if attr in cols:
            continue

    hdu = BinTableHDU(tbl, header=hdr, name=structure.__str__())

    return hdu
