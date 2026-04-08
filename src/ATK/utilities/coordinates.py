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

from ..configuration.survey_config import SURVEY_CONFIG
from ..queries.vizier.vizier_query import gaia_query_by_source
from ..structures.DataSet import DataSet
from ..structures.Target import Target
from ..utilities.defaults import RETURNS

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


def correct_coords(coords: list[SkyCoord], target_epoch: Time, get_correction: bool = False):
    """
    Not vectorised, but can't do this since astropy doesn't allow nan + non-nan distance/pm information in a non-scalar SkyCoord
    """

    if not isinstance(coords, list):
        coords = [coords]

    correction = []
    out_coords = []
    for c in coords:
        dt = (target_epoch - c.obstime).to_value(u.yr) * u.yr
        pm_valid = np.isfinite(c.pm_ra_cosdec.to_value(u.mas / u.yr)) and np.isfinite(c.pm_dec.to_value(u.mas / u.yr))
        dist_valid = c.distance.unit.physical_type == "length" and np.isfinite(c.distance.value)

        if not pm_valid:
            out_coords.append(c)
            correction.append("none")
            continue

        if dist_valid:
            new_c = c.apply_space_motion(target_epoch)
            correction.append("full")
        else:
            dec_rad = np.deg2rad(c.dec.to_value(u.deg))
            new_ra = c.ra.to(u.deg) + (c.pm_ra_cosdec.to(u.deg / u.yr) * dt / np.cos(dec_rad))
            new_dec = c.dec.to(u.deg) + (c.pm_dec.to(u.deg / u.yr) * dt)

            new_c = SkyCoord(ra=new_ra, dec=new_dec, pm_ra_cosdec=c.pm_ra_cosdec, pm_dec=c.pm_dec, frame=c.frame, obstime=target_epoch)

            correction.append("partial")

        out_coords.append(new_c)

        # if not np.isfinite(new_c.pm_ra_cosdec.value):
        #     print("IN:\n", c)
        #     print("OUT:\n", new_c, "\n\n\n")

    if get_correction:
        return out_coords, correction
    else:
        return out_coords


def get_gaia_target(source: int) -> Target:
    """
    Generates a SkyCoord using Gaia astrometry
    """

    gaia_epoch = SURVEY_CONFIG._get_epochs("vizier")["gaia"]

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

    return Target(copy.deepcopy(coord), copy.deepcopy(coord), None, source, "gaia", correction)


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
        epochs = SURVEY_CONFIG._get_epochs(query_kind)

        # If no epoch definition, can't correct
        if survey not in epochs:
            target.correction = "none"
            return target
        else:
            survey_epoch = epochs[survey]
    else:
        survey_epoch = epoch

    # correct coordinates to survey
    new_coords = correct_coords(target.coords, survey_epoch)[0]

    if not make_copy:
        target.coords = new_coords
        return target
    else:
        return replace(target, coords=new_coords)


def prepare_search(targets: list[Target], query_kind: str, survey: str = None, epoch: Time = None, **kwargs) -> tuple[SkyCoord, any]:
    """
    Prepares a search with an input skycoord. Returns the position of the search and a partially completed ATK DataSet
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
    structure = DataSet(kind=query_kind, targets=targets, exception=False)

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

    def correction(self):
        import math

        self.ra += (
            (self.year_delta * self.pmra / 3600000 + self.month_delta * self.pmra / 43200000) * 1 / math.cos(self.dec / 360 * 2 * math.pi)
        )
        self.dec += self.year_delta * self.pmdec / 3600000 + self.month_delta * self.pmdec / 43200000

        return [self.ra, self.dec]

    return expanded_radius
