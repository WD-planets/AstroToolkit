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


# ----
# MAIN
# ----


def write_structure(structure: any, path: str | Path):
    hdul = HDUList()

    query_hdu = struct_to_hdu(structure, ignore_attrs=("data", "frame", "epoch"), kind=PrimaryHDU)
    hdul.append(query_hdu)

    for attr, val in structure.__dict__.items():
        if attr == "data":
            if isinstance(val, pd.DataFrame):
                tbl = Table.from_pandas(val)
                hdu = BinTableHDU(tbl, header=Header(), name=structure.__repr__())
                hdul.append(hdu)
            elif isinstance(val, list):
                for ctr in val:
                    hdul.append(ctr.to_hdu())
            else:
                raise ValueError(f"Unexpected type of .data attribute in structure '{type(structure)}'.")

    hdul.writeto(path, overwrite=True)
