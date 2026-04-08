import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU, Header
from astropy.table import Table, hstack, vstack
from astropy.table.row import Row
from astropy.time import Time

from ..structures.DataSet import DataSet
from ..structures.Target import Target

# ============================================================
# SKYCOORD → TABLE
# ============================================================


def skycoord_to_table(coord: SkyCoord, prefix: str = ""):
    if coord.data.differentials:
        pmra = coord.pm_ra_cosdec.to_value(u.mas / u.yr)
        pmdec = coord.pm_dec.to_value(u.mas / u.yr)
    else:
        pmra = np.nan
        pmdec = np.nan

    if coord.distance != u.one:
        dist = coord.distance.to_value(u.pc)
    else:
        dist = np.nan

    row = {
        f"{prefix}ra": coord.ra.to_value(u.deg) * u.deg,
        f"{prefix}dec": coord.dec.to_value(u.deg) * u.deg,
        f"{prefix}frame": coord.frame.name,
        f"{prefix}epoch": coord.obstime.fits,
        f"{prefix}pm_ra_cosdec": pmra * u.mas / u.yr,
        f"{prefix}pm_dec": pmdec * u.mas / u.yr,
        f"{prefix}distance": dist * u.pc,
    }

    table = Table([row])

    return table


def targets_to_hdu(targets: list[Target]) -> BinTableHDU:
    init_coords = [skycoord_to_table(t.initial_coords, "input_") for t in targets]
    final_coords = [skycoord_to_table(t.coords, "final_") for t in targets]

    init_tbl = vstack(init_coords)
    final_tbl = vstack(final_coords)
    combined_tbl = hstack([init_tbl, final_tbl])

    combined_tbl["identifier"] = np.array([str(t.identifier) for t in targets])
    combined_tbl["survey"] = np.array([str(t.survey) for t in targets])
    combined_tbl["correction"] = np.array([t.correction for t in targets])
    combined_tbl["radius"] = [t.radius if t.radius is not None else np.nan for t in targets]

    return BinTableHDU(combined_tbl, header=Header(), name="TARGETING INFO")


def read_skycoord(row: Row, prefix: str = ""):
    ra = row[f"{prefix}ra"] * row.table[f"{prefix}ra"].unit
    dec = row[f"{prefix}dec"] * row.table[f"{prefix}dec"].unit
    frame = row[f"{prefix}frame"]
    epoch = Time(row[f"{prefix}epoch"], format="fits")

    pmra = row[f"{prefix}pm_ra_cosdec"]
    pmdec = row[f"{prefix}pm_dec"]
    dist = row[f"{prefix}distance"]

    kwargs = dict(ra=ra, dec=dec, frame=frame, obstime=epoch)

    if not np.ma.is_masked(pmra) and not np.ma.is_masked(pmdec):
        kwargs["pm_ra_cosdec"] = pmra * row.table[f"{prefix}pm_ra_cosdec"].unit
        kwargs["pm_dec"] = pmdec * row.table[f"{prefix}pm_dec"].unit

    if not np.ma.is_masked(dist):
        kwargs["distance"] = dist * row.table[f"{prefix}distance"].unit

    return SkyCoord(**kwargs)


def get_targets_from_hdu(structure: DataSet, primary_hdu: BinTableHDU, target_hdu: BinTableHDU):
    primary_header = primary_hdu.header

    if primary_header.get("ATK_FRAME") and hasattr(structure, "frame"):
        structure.frame = primary_header["ATK_FRAME"]

    if primary_header.get("ATK_EPOCH") and hasattr(structure, "epoch"):
        structure.epoch = primary_header["ATK_EPOCH"]

    table = Table.read(target_hdu)
    targets = []

    for row in table:
        init_coord = read_skycoord(row, prefix="input_")
        final_coord = read_skycoord(row, prefix="final_")

        identifier = row["identifier"] if row["identifier"] != "None" else None
        survey = row["survey"] if row["survey"] != "None" else None
        radius = row["radius"] * row.table["radius"].unit if not np.ma.is_masked(row["radius"]) else None

        targets.append(
            Target(
                initial_coords=init_coord,
                coords=final_coord,
                radius=radius,
                identifier=identifier,
                survey=survey,
                correction=row["correction"],
            )
        )

    structure.targets = targets

    return structure
