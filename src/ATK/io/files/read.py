import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits.hdu import BinTableHDU, ImageHDU
from astropy.table import Table

from ...structures.definitions import QueryResult
from ...utilities.mapping import build_structure_map

SKYCOORD_KEYS = ("ATK_RA", "ATK_DEC", "ATK_PMRA", "ATK_PMDEC", "ATK_DISTANCE", "ATK_FRAME", "ATK_EPOCH")

# -----------------
# PARSING FUNCTIONS
# -----------------


def parse_header_skycoord(path: str | Path, hdr: Header) -> SkyCoord:
    """
    Parses multiple header keys into an astropy SkyCoord object
    """

    def get_quantity(hdr, key, unit):
        val = hdr.get(key)
        if not val:
            return val
        return val * unit

    if any([key not in hdr.keys() for key in SKYCOORD_KEYS]):
        raise ValueError(f"Invalid SkyCoord decomposition in header of {path}")

    coord = SkyCoord(
        ra=get_quantity(hdr, "ATK_RA", u.deg),
        dec=get_quantity(hdr, "ATK_DEC", u.deg),
        pm_ra_cosdec=get_quantity(hdr, "ATK_PMRA", u.mas / u.yr),
        pm_dec=get_quantity(hdr, "ATK_PMDEC", u.mas / u.yr),
        distance=get_quantity(hdr, "ATK_DISTANCE", u.pc),
        frame=hdr.get("ATK_FRAME"),
        obstime=hdr.get("ATK_EPOCH"),
    )

    return coord


def parse_primary_header(path: str | Path, hdr: Header) -> QueryResult:
    structure = QueryResult()

    if "ATK_SKYCOORD" not in hdr:
        raise ValueError(f"No targeting information found in primary header of {path}.")

    structure.position = parse_header_skycoord(path, hdr)

    structure = parse_header(path, hdr, structure)

    return structure


def parse_header(path: str | Path, hdr: Header, obj: object) -> object:
    ATK_keys = {key: val for key, val in hdr.items() if key.startswith("ATK_")}
    for key, val in ATK_keys.items():
        attr = key[4:].lower()

        if hasattr(obj, attr):
            setattr(obj, attr, val)

    return obj


def parse_generic_bintable(structure: QueryResult, path: str | Path, hdr: Header, data: any):
    structure_map = build_structure_map()

    df = Table(data).to_pandas()

    ctnr = structure_map.get(hdr.get("ATK_KIND"), lambda: None)()
    if not ctnr:
        structure.data = df
        return structure

    ctnr = parse_header(path, hdr, ctnr)

    for col in df:
        if hasattr(ctnr, col):
            setattr(ctnr, col, np.asarray(df[col]))

    structure.data = [ctnr]

    return structure


# ----
# MAIN
# ----


def read_local(path: str | Path):
    hdul = fits.open(path)

    structure = parse_primary_header(path, hdul[0].header)

    for hdu in hdul[1:]:
        if not hdu.header.get("ATK_EXT", None):
            warnings.warn(f"ATK: Non-ATK extension found in target file {path} has been ignored.")
            continue
        if not isinstance(hdu, (BinTableHDU, ImageHDU)):
            warnings.warn(f"ATK: Unexpected table type '{type(hdu)}' in target file {path} has been ignored.")
            continue

        structure = parse_generic_bintable(structure, path, hdu.header, hdu.data)

    structure.show()

    hdul.close()
