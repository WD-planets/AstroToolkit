import warnings
from pathlib import Path

import pandas as pd
from astropy.io.fits import HDUList, Header
from astropy.io.fits.hdu import BinTableHDU, PrimaryHDU
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table

from ...io.structure_io import struct_to_hdu
from ...utilities.misc import get_package_version
from ..target_io import targets_to_hdu

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

    # .kind handled in structure_io, targets handled below
    query_hdu = struct_to_hdu(structure, ignore_attrs=["kind", "targets", "data", "figure"], hdu_kind=PrimaryHDU)
    hdul.append(query_hdu)

    if hasattr(structure, "targets"):
        hdul.append(targets_to_hdu(structure.targets))

    # iterate through .data
    for ctr in structure.data:
        # write dataframe to hdu (e.g. in Vizier queries)
        if isinstance(ctr, pd.DataFrame):
            hdul.append(dataframe_to_hdu(structure, ctr))
            continue

        # otherwise use .to_hdu() method of ATK container
        hdus = ctr.to_hdu()
        if not isinstance(hdus, (tuple, list)):
            hdul.append(hdus)
            continue

        # if multiple hdus returned (e.g. images)
        for hdu in hdus:
            hdul.append(hdu)

    # store ATK version used to generate file
    primary_hdr = hdul[0].header
    primary_hdr.append(("ATK_VER", get_package_version(), "ATK version at time of file creation"))

    # tag all HDUs as coming from ATK
    for hdu in hdul:
        hdr = hdu.header
        hdr.append(("ATK_EXT", True, "If True, this is a fits file from ATK"))

    hdul.writeto(path, overwrite=True)

    return path
