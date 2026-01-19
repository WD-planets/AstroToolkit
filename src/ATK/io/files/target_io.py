from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU, Header
from astropy.table import Table, hstack
from astropy.table.row import Row
from astropy.time import Time
from astropy.units import Unit

from ...structures.definitions import BaseQueryResult, Target


def stack_skycoords(coords: list[SkyCoord]):
    """
    Combines multiple SkyCoords into a single SkyCoord while retaining all proper motion and frame information
    """

    frame = coords[0].frame.name
    obstime = coords[0].obstime

    return SkyCoord(
        ra=[c.ra for c in coords],
        dec=[c.dec for c in coords],
        pm_ra_cosdec=[c.pm_ra_cosdec for c in coords],
        pm_dec=[c.pm_dec for c in coords],
        distance=[c.distance for c in coords],
        frame=frame,
        obstime=obstime,
    )


def skycoord_to_table(coord: SkyCoord, prefix: str = ""):
    table = Table()

    table[f"{prefix}ra"] = coord.ra
    table[f"{prefix}dec"] = coord.dec
    table[f"{prefix}pm_ra_cosdec"] = coord.pm_ra_cosdec
    table[f"{prefix}pm_dec"] = coord.pm_dec
    table[f"{prefix}distance"] = coord.distance
    table[f"{prefix}frame"] = coord.frame.name
    table[f"{prefix}epoch"] = coord.obstime.fits

    return table


def targets_to_hdu(targets: list[Target]) -> BinTableHDU:
    """
    Converts a list of targets to a BinTableHDU
    """

    # stack all initial Target coords and convert to dataframe
    init_coords = stack_skycoords([target.initial_coords for target in targets])
    init_coords_tbl = skycoord_to_table(init_coords, prefix="input_")

    # stack all final Target coords and convert to dataframe
    final_coords = stack_skycoords([target.coords for target in targets])
    final_coords_tbl = skycoord_to_table(final_coords, prefix="final_")

    combined_tbl = hstack([init_coords_tbl, final_coords_tbl])

    # add basic info
    combined_tbl["identifier"] = [target.identifier for target in targets]
    combined_tbl["survey"] = [target.survey for target in targets]
    combined_tbl["correction"] = [target.correction for target in targets]

    combined_tbl.remove_columns(["final_pm_ra_cosdec", "final_pm_dec", "final_distance"])
    combined_tbl.rename_columns(
        ["input_pm_ra_cosdec", "input_pm_dec", "input_distance"], ["pm_ra_cosdec", "pm_dec", "distance"]
    )

    hdu = BinTableHDU(combined_tbl, header=Header(), name="TARGETING INFO")

    return hdu


def get_unit(row: Row, column: str) -> Unit:
    return row.table[column].unit


def read_skycoord(row: Row, prefix: str = ""):
    ra = row[f"{prefix}ra"] * get_unit(row, f"{prefix}ra")
    dec = row[f"{prefix}dec"] * get_unit(row, f"{prefix}dec")
    pmra = row["pm_ra_cosdec"] * get_unit(row, "pm_ra_cosdec")
    pmdec = row["pm_dec"] * get_unit(row, "pm_dec")
    distance = row["distance"] * get_unit(row, "distance")
    frame = row[f"{prefix}frame"]
    epoch = row[f"{prefix}epoch"]

    coord = SkyCoord(
        ra=ra,
        dec=dec,
        pm_ra_cosdec=pmra,
        pm_dec=pmdec,
        distance=distance,
        frame=frame,
        obstime=Time(epoch, format="fits"),
    )

    return coord


def get_targets_from_hdu(structure: BaseQueryResult, primary_hdu: BinTableHDU, target_hdu: BinTableHDU) -> list[Target]:
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

        target = Target(init_coord, final_coord, row["identifier"], row["survey"], row["correction"])
        targets.append(target)

    structure.targets = targets

    return structure
