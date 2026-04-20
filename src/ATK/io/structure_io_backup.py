from __future__ import annotations

import typing
from dataclasses import fields, is_dataclass
from types import NoneType, UnionType
from typing import TYPE_CHECKING, get_args, get_origin

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU, PrimaryHDU
from astropy.table import Table
from astropy.units import Quantity

from ..structures.Target import Target

if TYPE_CHECKING:
    from ..structures.Image import Image
    from ..structures.structures_core import Container

# types (in typehints) that should be considered as being columns of a dataframe
COLUMN_TYPES = (np.ndarray, pd.Series)

BASIC_TYPES = (int, float, str, bool, NoneType)


def is_scalar_quantity(val):
    return isinstance(val, Quantity) and getattr(val, "shape", None) == ()


# ----------------------
# HEADER WRITE FUNCTIONS
# ----------------------


def write_fallback(attr: str, hdr: Header, key: str, value: any) -> Header:
    """
    Writes generic data types to the header (e.g. str, int, float)
    """

    try:
        hdr.append((f"ATK_{key.upper()}", value))
        return hdr
    except Exception:
        pass

    try:
        hdr.append((f"ATK_{key.upper()}", str(value)))
        return hdr
    except Exception:
        raise ValueError(
            f"Failed to write value '{value}' of type '{type(value)}' in attribute '{attr}' to FITS header key 'ATK_{key.upper()}'."
        )

    return hdr


def write_quantity(attr: str, hdr: Header, key: str, value: Quantity) -> Header:
    """
    Writes an astropy quantity to a fits header as two keys, ATK_... with unit ATK_..._U
    """

    write_fallback(attr, hdr, key, value.value)

    hdr.append((f"ATK_{key.upper()}_U", value.unit.to_string("fits")))

    return hdr


def write_skycoord(attr: str, hdr: Header, coord: SkyCoord) -> Header:
    """
    Splits a SkyCoord into its components and writes this as a set of header keys
    """

    hdr.append(("ATK_SC", attr, "Origin attribute of decomposed SkyCoord"))
    hdr.append(("ATK_FRAME", coord.frame.name, "Coordinate frame of ATK_RA and ATK_DEC"))

    if coord.obstime:
        hdr.append(("ATK_EPOCH", coord.obstime.fits, "Epoch of ATK_RA and ATK_DEC"))
    else:
        hdr.append(("ATK_EPOCH", None))

    # proper motion data
    if coord.data.differentials:
        hdr = write_quantity(attr, hdr, "pmra", coord.pm_ra_cosdec)
        hdr = write_quantity(attr, hdr, "pmdec", coord.pm_dec)
    else:
        hdr.append(("ATK_PMRA", None))
        hdr.append(("ATK_PMDEC", None))

    # distance
    if coord.distance != u.one:
        hdr = write_quantity(attr, hdr, "distance", coord.distance)
    else:
        hdr.append(("ATK_DISTANCE", None))

    hdr = write_quantity(attr, hdr, "ra", coord.ra)
    hdr = write_quantity(attr, hdr, "dec", coord.dec)

    return hdr


# Map to special header writing functions
WRITE_MAP = {
    SkyCoord: lambda **kwargs: write_skycoord(kwargs["attr"], kwargs["hdr"], kwargs["value"]),
    Quantity: lambda **kwargs: write_quantity(kwargs["attr"], kwargs["hdr"], kwargs["key"], kwargs["value"]),
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
    """
    Returns a list of array-like attributes of a data structure using typehinting
    """

    # ensure we have the class type, not an instance
    cls = structure if isinstance(structure, type) else structure.__class__
    if not is_dataclass(cls):
        raise TypeError("get_cols expects a dataclass.")

    cols = []
    hints = typing.get_type_hints(cls)

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


def struct_to_table(structure: any) -> Table:
    """
    Combines array-like attributes of a data structure into a single astropy Table
    """

    cols = get_cols(structure)

    table = Table()
    for col in cols:
        val = getattr(structure, col)
        if val is None:
            continue

        if not isinstance(val, COLUMN_TYPES):
            val = [val]

        if isinstance(val, Quantity) and is_scalar_quantity(val):
            table[col] = val
        else:
            table[col] = val

    return table


def struct_from_table(ctnr: any, target: Target | int | SkyCoord, data: Table, **kwargs: dict) -> any:
    from ..queries.query_core import setup_targeting

    ctnr_cols = get_cols(ctnr)

    relevant_data = {}
    for col_name in data.colnames:
        if hasattr(ctnr, col_name) and col_name in ctnr_cols:
            col = data[col_name]
            unit = getattr(col, "unit")

            if unit is not None:
                relevant_data[col_name] = Quantity(col, unit=unit)
            else:
                relevant_data[col_name] = np.array(col)

    for arg, val in kwargs.items():
        if hasattr(ctnr, arg) and arg not in ctnr_cols:
            relevant_data[arg] = val

    missing = [f"'{param}'" for param in ctnr._required if param not in kwargs]
    if missing:
        raise ValueError(f"Missing required arguments(s) for {ctnr.__name__} construction: {', '.join(missing)}")

    ctnr = ctnr(**relevant_data)

    targets = setup_targeting(target)
    if len(targets) > 1:
        raise ValueError("Only one target may be provided per data container.")

    ctnr._target_key = targets[0]._key

    return ctnr


def struct_to_dataframe(structure: any) -> pd.DataFrame:
    """
    Combines the array-like attributes of a data structure into a single pandas DataFrame
    """

    cols = get_cols(structure)

    data = {}
    for col in cols:
        val = getattr(structure, col)
        if val is None:
            continue

        if not isinstance(val, COLUMN_TYPES):
            val = [val]

        print(col, val)
        if isinstance(val, Quantity) and val.isscalar:
            print("skipped")
            continue

        data[col] = val

    return pd.DataFrame.from_dict(data)


def struct_from_dataframe(ctnr: any, target: Target | int | SkyCoord, data: pd.DataFrame, **kwargs) -> any:
    from ..queries.query_core import setup_targeting

    ctnr_cols = get_cols(ctnr)

    relevant_data = {}
    for col in data.columns.values.tolist():
        if hasattr(ctnr, col) and col in ctnr_cols:
            relevant_data[col] = data[col].to_numpy()

    for arg, val in kwargs.items():
        if hasattr(ctnr, arg) and arg not in ctnr_cols:
            relevant_data[arg] = val

    missing = [f"'{param}'" for param in ctnr._required if param not in kwargs]
    if missing:
        raise ValueError(f"Missing required arguments(s) for {ctnr.__name__} construction: {', '.join(missing)}")

    ctnr = ctnr(**relevant_data)

    targets = setup_targeting(target)
    if len(targets) > 1:
        raise ValueError("Only one target may be provided per data container.")

    ctnr._target_key = targets[0]._key

    return ctnr


def struct_to_hdu(structure: any, ignore_attrs: list = [], hdu_kind: BinTableHDU | PrimaryHDU = BinTableHDU) -> BinTableHDU:
    """
    Convert a data structure into a fits HDU.
    """

    hdr = Header()

    # combine array-like attributes into a dataframe
    cols = get_cols(structure)
    tbl = struct_to_table(structure)

    # PrimaryHDU stores query kind, extensions store data container kind
    kind = structure.__dict__.get("kind", type(structure).__name__)
    hdr.append(("ATK_KIND", kind))

    # write target key for mapping targets to containers (exception to ignoring attrs with _ below)
    if hasattr(structure, "_target_key"):
        hdr.append(("ATK__TARGET_KEY", structure._target_key))

    # ^
    if hasattr(structure, "_cache_key"):
        hdr.append(("ATK__CACHE_KEY", structure._cache_key))

    # iterate through remaining structure attributes and write them to the header
    for attr, val in structure.__dict__.items():
        if attr.startswith("_"):
            continue
        if attr in ignore_attrs:
            continue
        if attr in cols:
            continue

        if isinstance(val, Quantity):
            if val.isscalar:
                hdr = write_quantity(attr, hdr, attr, val)
            else:
                tbl[attr] = val
        elif type(val) in WRITE_MAP:
            hdr = WRITE_MAP[type(val)](attr=attr, hdr=hdr, key=attr, value=val)
        else:
            hdr = write_fallback(attr, hdr, attr, val)

    # generate HDU
    extname = structure.__str__().lstrip("<").rstrip(">")

    tbl = normalise_table_strings(tbl)

    if hdu_kind == BinTableHDU:
        hdu = BinTableHDU(tbl, header=hdr, name=extname)
    elif hdu_kind == PrimaryHDU:
        hdu = PrimaryHDU(tbl, header=hdr)
    else:
        raise ValueError(f"Unexpected fits extension '{hdu_kind}'.")

    return hdu


# -----------------------
# SPECIAL TRANSFORMATIONS
# -----------------------


def normalise_table_strings(table: Table):
    for name in table.colnames:
        col = table[name]

        if col.dtype.kind in {"U", "S", "O"}:
            arr = np.array(col, dtype=str)

            # replace masked/null-like values
            arr[arr == "None"] = ""
            arr[arr == "masked"] = ""

            maxlen = max(len(x) for x in arr) if len(arr) > 0 else 1
            table[name] = arr.astype(f"<U{maxlen}")

    return table


def image_to_hdu(image: Image):
    """
    Converts an ATK image to a HDU, needed as an Image is itself a hdu
    """

    hdu = image.hdu
    hdr = hdu.header

    hdr.append(("ATK_KIND", "Image", "ATK container kind"))
    for attr, val in image.__dict__.items():
        # basic types and those with special writing functions (e.g. SkyCoord decomposition)
        if type(val) in WRITE_MAP or type(val) in BASIC_TYPES:
            hdr = WRITE_MAP.get(type(val), write_fallback)(attr=attr, hdr=hdr, key=attr, value=val)

    extname = image.__str__().lstrip("<").rstrip(">")
    image_hdu = ImageHDU(data=hdu.data, header=hdr, name=extname)

    overlay_hdr = Header()
    overlay_hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))
    overlay = image.overlay
    if overlay is None:
        overlay = pd.DataFrame()
    table = Table.from_pandas(overlay)
    overlay_hdu = BinTableHDU(data=table, header=overlay_hdr, name="<Overlay Data>")

    return (image_hdu, overlay_hdu)


def simple_to_hdu(entry: Container):
    hdr = Header()

    hdr.append(("ATK_KIND", type(entry).__name__, "ATK container kind"))
    hdr.append(("ATK_SIMPLE", True, "If True, data originates from as a table"))

    for attr, val in entry.__dict__.items():
        if type(val) in WRITE_MAP or type(val) in BASIC_TYPES:
            hdr = WRITE_MAP.get(type(val), write_fallback)(attr=attr, hdr=hdr, key=attr, value=val)

    extname = entry.__str__().lstrip("<").rstrip(">")

    # normalise string lengths
    table = entry.table
    table = normalise_table_strings(table)

    hdu = BinTableHDU(data=entry.table, header=hdr, name=extname)

    return hdu
