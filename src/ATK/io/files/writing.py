from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import HDUList, Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.table import Table
from astropy.time import Time

WRITE_MAP = {}


# ---------------
# WRITE FUNCTIONS
# ---------------


# -------
# MAPPING
# -------


WRITE_MAP.update(
    {
        pd.DataFrame: write_dataframe,
        SkyCoord: lambda data, name, hdr, hdul: write_skycoord(data, hdr, hdul),
        dict: lambda data, name, hdul, hdr: write_dataframe(data.to_dict(orient="list"), name, hdul, hdr),
    }
)


# ----
# MAIN
# ----


def write_structure(structure: any, fname: str | Path):
    hdr = Header()
    hdul = HDUList()

    for attr, val in structure.__dict__.items():
        write_function = WRITE_MAP.get(type(val), None)
        if not write_function:
            raise Exception(f"No write function found for dtype '{type(val)}'.")

        hdr, hdul = write_function()
