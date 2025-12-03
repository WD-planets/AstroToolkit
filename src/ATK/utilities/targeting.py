import warnings

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from erfa import ErfaWarning

from ..configuration.epoch_config import EPOCH_CONFIG
from ..queries.vizier.vizier_query import gaia_query_by_source
from ..utilities.defaults import RETURNS
from ..utilities.mapping import build_structure_map

# ignore bad distance warning
warnings.filterwarnings("ignore", category=ErfaWarning)


def check_distance(coord: SkyCoord) -> bool:
    """
    Checks if a SkyCoord has a distance
    """

    return coord.distance.unit.is_equivalent(u.pc)


def check_finite(val: any) -> bool:
    """
    Wrapper for np.isfinite to also include nans
    """

    return np.isfinite(val) if val is not None else False


def check_targeting(target: int | SkyCoord) -> tuple[int | None, SkyCoord | None]:
    """
    Checks that a position or source ID have been provided, and sets up a SkyCoord with the expected format in the former case
    """

    if isinstance(target, int):
        source = target
        position = None
    elif isinstance(target, SkyCoord):
        source = None
        position = target

    if position:
        if position.frame != "icrs":
            icrs_position = position.transform_to("icrs")
        else:
            icrs_position = position
        icrs_position.obstime
        icrs_ra = icrs_position.ra
        icrs_dec = icrs_position.dec
        time = Time("2000-01-01", format="iso")
        final_position = SkyCoord(icrs_ra, icrs_dec, frame="icrs", obstime=time)
    else:
        final_position = None

    return source, final_position


def get_gaia_skycoord(source: int) -> SkyCoord:
    """
    Converts Gaia data into a SkyCoord containing all necessary astrometry
    """

    gaia_epoch = EPOCH_CONFIG.as_dict()["vizier_aliases"]["gaia"]

    gaia_data = gaia_query_by_source(source)
    if gaia_data is RETURNS.EXCEPTION or gaia_data is RETURNS.NULL:
        return gaia_data

    ra, dec = gaia_data["RA_ICRS"].tolist()[0], gaia_data["DE_ICRS"].tolist()[0]
    pmra, pmdec = gaia_data["pmRA"].tolist()[0], gaia_data["pmDE"].tolist()[0]
    parallax = gaia_data["Plx"].tolist()[0]

    if not np.isfinite(ra) or not np.isfinite(dec):
        raise Exception(f"Couldn't get RA and DEC for Gaia source '{source}'.")

    # if pmra or pmdec are bad, don't supply these to SkyCoord -> no correction
    if not check_finite(pmra) or not check_finite(pmdec):
        pmra, pmdec = None, None
    else:
        pmra, pmdec = pmra * u.mas / u.yr, pmdec * u.mas / u.yr

    # if parallax is bad or negative, do not supply a distance to SkyCoord -> partial correction
    if not check_finite(parallax) or parallax <= 0:
        distance = None
    else:
        distance = 1000 / parallax * u.pc

    coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, pm_ra_cosdec=pmra, pm_dec=pmdec, distance=distance, obstime=gaia_epoch, frame="icrs")

    return coord


def correct_skycoord(
    position: SkyCoord, query_kind: str = None, survey: str = None, epoch: Time = None, get_correction_degree: bool = False
) -> SkyCoord:
    """
    Corrects a SkyCoord to a given epoch definition from a given section or a given epoch, returns corrected SkyCoord object and optionally returns the success level of the correction. If correction isn't possible, just returns the original SkyCoord
    """

    if (survey is None) == (epoch is None):
        raise ValueError("Specify exactly one of 'survey', 'epoch'.")

    correction_degree = "full"

    if survey:
        epochs = EPOCH_CONFIG.get_section_by_query_kind(query_kind)

        # If no epoch definition, can't correct
        if survey not in epochs:
            corrected_position = position
            correction_degree = "none"
        else:
            survey_epoch = epochs[survey]
    else:
        survey_epoch = epoch

    # if proper motion data missing, can't correct
    if not position.data.differentials and correction_degree != "none":
        corrected_position = position
        correction_degree = "none"

    # if distance missing due to bad parallax, partial correction (usually fine)
    if position.distance == u.one and correction_degree != "none":
        correction_degree = "partial"

    # correct coordinates to survey
    if correction_degree != "none":
        corrected_position = position.apply_space_motion(survey_epoch)

    if get_correction_degree:
        return corrected_position, correction_degree
    else:
        return corrected_position


def prepare_search(target: int | SkyCoord, query_kind: str, survey: str = None, epoch: Time = None, **kwargs) -> tuple[SkyCoord, any]:
    """
    Prepares a search with an input position or source. Returns the position of the search and a partially completed data structure
    """

    source, position = check_targeting(target)

    if source:
        # translate Gaia position to epoch of survey if an epoch definition exists
        gaia_pos = get_gaia_skycoord(source)
        if gaia_pos in [RETURNS.NULL, RETURNS.EXCEPTION]:
            search_pos = None
            correction_degree = "none"
        elif survey != "gaia":
            search_pos, correction_degree = correct_skycoord(gaia_pos, query_kind, survey=survey, epoch=epoch, get_correction_degree=True)
        else:
            search_pos = gaia_pos
            correction_degree = "n/a"

    elif position:
        # needed below when checking if an exception occured at this point
        gaia_pos = None

        correction_degree = "none"
        search_pos = position

    structure_map = build_structure_map()
    structure = structure_map[query_kind](
        kind=query_kind,
        survey=survey,
        position=search_pos,
        source=source,
        radius=kwargs.get("radius", None),
        frame=getattr(getattr(search_pos, "frame", None), "name", None),
        epoch=getattr(search_pos, "obstime", None),
        correction=correction_degree,
        exception=True if gaia_pos == RETURNS.EXCEPTION else False,
    )

    return search_pos, structure


def correct_radius(target: int | SkyCoord, radius: float, epoch: Time = None):
    source, position = check_targeting(target)

    # return uncorrected radius
    if position:
        return radius

    gaia_pos = get_gaia_skycoord(source)
    corrected_pos = correct_skycoord(gaia_pos, epoch=epoch)

    separation = corrected_pos.separation(gaia_pos)
    expanded_radius = (radius * u.arcsec + separation.to(u.arcsec)).value

    return expanded_radius
