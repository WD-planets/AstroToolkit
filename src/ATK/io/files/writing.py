import warnings
from pathlib import Path

import pandas as pd
from astropy.io.fits import HDUList, Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table

from ...structures.structure_io import BASIC_TYPES, struct_to_hdu

warnings.simplefilter("ignore", category=VerifyWarning)

WRITE_MAP = {}

# -----------------------
# ADDITIONAL TRANSLATIONS
# -----------------------


def dataframe_to_hdu(structure, data: pd.DataFrame) -> BinTableHDU:
    """
    Convert a pandas dataframe to an astropy BinTableHDU
    """

    hdr = Header()
    hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))
    tbl = Table.from_pandas(data)
    hdu = BinTableHDU(tbl, header=hdr, name=structure.__repr__())

    return hdu


# ----
# MAIN
# ----


def write_local(structure: any, path: str | Path) -> Path:
    """
    Write an ATK data structure to a local fits file
    """

    path = path or structure._fname

    hdul = HDUList()

    query_hdu = struct_to_hdu(structure, ignore_attrs=["data", "frame", "epoch"], hdu_kind=PrimaryHDU)
    hdul.append(query_hdu)

    # iterate through structure attributes
    for attr, val in structure.__dict__.items():
        if isinstance(val, BASIC_TYPES):
            continue

        if not isinstance(val, list):
            continue

        for ctr in val:
            if isinstance(ctr, pd.DataFrame):
                hdul.append(dataframe_to_hdu(structure, ctr))
                continue

            hdus = ctr.to_hdu()
            if not isinstance(hdus, (tuple, list)):
                hdul.append(hdus)
                continue

            for hdu in hdus:
                hdul.append(hdu)

    hdul.writeto(path, overwrite=True)

    return path
