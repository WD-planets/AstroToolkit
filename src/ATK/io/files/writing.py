import warnings
from pathlib import Path

import pandas as pd
from astropy.io.fits import HDUList, Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table

from ...structures.structure_io import struct_to_hdu

warnings.simplefilter("ignore", category=VerifyWarning)

WRITE_MAP = {}

# -----------------------
# ADDITIONAL TRANSLATIONS
# -----------------------


def dataframe_to_hdu(structure, data: pd.DataFrame) -> BinTableHDU:
    hdr = Header()
    hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))
    tbl = Table.from_pandas(data)
    hdu = BinTableHDU(tbl, header=hdr, name=structure.__repr__())

    return hdu


# ----
# MAIN
# ----


def write_local(structure: any, path: str | Path) -> Path:
    path = path or structure._fname

    hdul = HDUList()

    query_hdu = struct_to_hdu(structure, ignore_attrs=["data", "frame", "epoch"], hdu_kind=PrimaryHDU)
    hdul.append(query_hdu)

    for attr, val in structure.__dict__.items():
        if attr == "data":
            if isinstance(val, pd.DataFrame):
                hdul.append(dataframe_to_hdu(structure, val))
            elif isinstance(val, list):
                for ctr in val:
                    hdul.append(ctr.to_hdu())
            else:
                raise ValueError(f"Unexpected type of .data attribute in structure '{type(structure)}'.")

    hdul.writeto(path, overwrite=True)

    return path
