import copy
import warnings
from dataclasses import replace

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from erfa import ErfaWarning

from ..configuration.epoch_config import EPOCH_CONFIG
from ..queries.vizier.vizier_query import gaia_query_by_source
from ..structures.definitions import Target
from ..utilities.defaults import RETURNS
from ..utilities.mapping import get_query_result_map

# ignore bad distance warning
warnings.filterwarnings("ignore", category=ErfaWarning)

REQUIRED_COLS = ["ra", "dec", "pm_ra_cosdec", "pm_dec"]


def check_finite(val: any) -> bool:
    """
    Wrapper for np.isfinite to also include None
    """

    return np.isfinite(val) if val is not None else False


def check_correction(coord: SkyCoord) -> str:
    """
    Gets the maximum possible correction degree for an ATK target based on its astrometry
    """

    if not coord.data.differentials:
        return "none"
    if coord.distance == u.one:
        return "partial"
    else:
        return "full"


def get_gaia_target(source: int) -> Target:
    """
    Generates a SkyCoord using Gaia astrometry
    """

    gaia_epoch = EPOCH_CONFIG.as_dict()["vizier_aliases"]["gaia"]

    gaia_data = gaia_query_by_source(source)
    if gaia_data is RETURNS.EXCEPTION or gaia_data is RETURNS.NULL:
        return gaia_data

    ra, dec = gaia_data["RA_ICRS"].tolist()[0], gaia_data["DE_ICRS"].tolist()[0]
    pmra, pmdec = gaia_data["pmRA"].tolist()[0], gaia_data["pmDE"].tolist()[0]
    parallax = gaia_data["Plx"].tolist()[0]

    if not np.isfinite(ra) or not np.isfinite(dec):
        return RETURNS.NULL

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

    correction = check_correction(coord)

    return Target(copy.deepcopy(coord), copy.deepcopy(coord), source, "gaia", correction)


def correct_target(target: Target, survey: str = None, epoch: Time = None, query_kind: str = None, make_copy=False) -> Target:
    """
    Corrects the SkyCoord of an ATK Target to a given epoch definition from a given section or a given epoch, returns corrected Target (if correction was possible)
    If make_copy is True, corrected Target is returned as a new instance (leaving the original unchanged)
    """

    if (survey is None) == (epoch is None):
        raise ValueError("Specify exactly one of 'survey', 'epoch'.")

    if target.correction == "none":
        return target

    if survey:
        epochs = EPOCH_CONFIG.get_section_by_query_kind(query_kind)

        # If no epoch definition, can't correct
        if survey not in epochs:
            return target
        else:
            survey_epoch = epochs[survey]
    else:
        survey_epoch = epoch

    # correct coordinates to survey
    new_coords = target.coords.apply_space_motion(survey_epoch)

    if not make_copy:
        target.coords = new_coords
        return target
    else:
        return replace(target, coords=new_coords)


def prepare_search(targets: list[Target], query_kind: str, survey: str = None, epoch: Time = None, **kwargs) -> tuple[SkyCoord, any]:
    """
    Prepares a search with an input skycoord. Returns the position of the search and a partially completed ATK QueryResult
    """

    corrected_targets = []
    for target in targets:
        # if a catalogue is provided for Vizier queries, treat this as a survey for correction
        if kwargs.get("catalogue"):
            survey = kwargs["catalogue"]

        # correct target to given survey/epoch
        if not kwargs.get("defer_correction", False):
            corrected_target = correct_target(target, survey, epoch, query_kind)
        else:
            corrected_target = target
        corrected_targets.append(corrected_target)

    # create requested structure
    structure_map = get_query_result_map()
    structure = structure_map[query_kind](
        kind=query_kind,
        survey=survey,
        targets=targets,
        radius=kwargs.get("radius", None),
        frame=getattr(
            getattr(corrected_targets[0].initial_coords, "frame", None), "name", None
        ),  # corrected targets should all have same frame
        epoch=getattr(corrected_targets[0].initial_coords, "obstime", None),  # corrected targets should all have same epoch
        correction=target.correction,
        exception=False,
    )

    return corrected_targets, structure


def correct_radius(target: Target, radius: Quantity, query_kind: str, survey: str):
    """
    Expands a search radius for proper motion
    """

    # return uncorrected radius if target can't be corrected
    if target.correction == "none":
        return radius

    corrected_target = correct_target(target, survey=survey, query_kind=query_kind, make_copy=True)
    separation = corrected_target.coords.separation(target.coords)
    expanded_radius = radius + separation.to(radius.unit)

    return expanded_radius


def dataframe_to_skycoord(data: pd.DataFrame, epoch: Time):
    """
    Converts a DataFrame containing any 'ra', 'dec', 'pm_ra_cosdec', and 'pm_dec' columns to a single SkyCoord
    """

    coords = SkyCoord(
        ra=data["ra"].to_numpy() * u.deg,
        dec=data["dec"].to_numpy() * u.deg,
        pm_ra_cosdec=data["pm_ra_cosdec"].to_numpy() * (u.mas / u.yr),
        pm_dec=data["pm_dec"].to_numpy() * (u.mas / u.yr),
        frame="icrs",
        obstime=epoch,
    )

    return coords


def skycoord_to_dataframe(coord: SkyCoord):
    """
    Converts a SkyCoord containing any number of positions to a DataFrame with 'ra', 'dec', 'pm_ra_cosdec' and 'pm_dec' columns
    """

    df = pd.DataFrame()

    df["ra"] = coord.ra.deg
    df["dec"] = coord.dec.deg
    df["pm_ra_cosdec"] = coord.pm_ra_cosdec.to(u.mas / u.yr).value
    df["pm_dec"] = coord.pm_dec.to(u.mas / u.yr).value

    return df


def correct_skycoord(coord: SkyCoord, input_epoch: Time, target_epoch: Time):
    """
    Corrects the positions in a SkyCoord for proprer motion
    """

    df = skycoord_to_dataframe(coord)
    df = correct_dataframe_coords(df, input_epoch, target_epoch)
    corrected_coord = dataframe_to_skycoord(df, input_epoch)

    return corrected_coord


def correct_dataframe_coords(data: pd.DataFrame, input_epoch: Time, target_epoch: Time, output_cols: list = []):
    """
    Corrects the coordinates in the 'ra' and 'dec' columns of a dataframe for proper motion in corresponding 'pm_ra_cosdec' and 'pm_dec' columns.
    Optionally saves the resulting coordinates to two new output columns (output_cols)
    """

    for col in REQUIRED_COLS:
        if col not in data:
            raise ValueError(f"DataFrame missing required column '{col}'.")

    # get mask of values with invalid PM information
    bad_pm_mask = data["pm_ra_cosdec"].isna() & data["pm_dec"].isna()

    # create output ra and dec columns if needed
    if output_cols:
        data[output_cols[0]] = data["ra"]
        data[output_cols[1]] = data["dec"]
    else:
        output_cols = ["ra", "dec"]

    good_pm_data = data.loc[~bad_pm_mask]
    if not good_pm_data.empty:
        # convert dataframe coordinates to skycoord
        coord = dataframe_to_skycoord(good_pm_data, input_epoch)

        # apply correction
        coord = coord.apply_space_motion(target_epoch)

        # update output ra and dec columns (only in rows where a correction occurred)
        data.loc[good_pm_data.index, output_cols[0]] = coord.ra.deg
        data.loc[good_pm_data.index, output_cols[1]] = coord.dec.deg

    return data
