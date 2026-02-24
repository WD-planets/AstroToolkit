import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.time import Time
from astropy.wcs.utils import proj_plane_pixel_scales

from ...configuration.base_config import BASE_CONFIG
from ...configuration.epoch_config import EPOCH_CONFIG
from ...configuration.overlay_config import OVERLAY_CONFIG
from ...structures.Image import Image
from ...structures.Target import Target
from ...Tools.query import query
from ...utilities.coordinates import correct_radius, correct_skycoord, dataframe_to_skycoord
from ...utilities.defaults import RETURNS
from ..simbad.simbad_query import get_ids


def get_overlay_data(image: Image, target: int | SkyCoord, survey: str, survey_info: dict, disable_corrections=False) -> pd.DataFrame:
    """
    Searches an image for any survey detections, and performs proper motion correction via piggybacking with astrometry from the chosen astrometric backend
    """

    # correct search radius (i.e. size of image) for maximum possible proper motion of object
    # between image epoch and non-gaia survey epoch. Padded by 25% to account for error
    radius = correct_radius(target, image.size, "vizier", survey) * 1.25
    piggyback_radius = BASE_CONFIG._get("overlay_settings", "piggyback_radius")

    # get non-gaia data
    non_gaia_data = query(kind="vizier", targets=image.search_pos, radius=radius, survey=survey).data
    if not non_gaia_data:
        return pd.DataFrame()
    else:
        non_gaia_data = non_gaia_data[0].data

    # get gaia data
    gaia_data = query(kind="vizier", targets=image.search_pos, radius=radius, survey="gaia").data
    if not gaia_data:
        gaia_data = pd.DataFrame()
    else:
        gaia_data = gaia_data[0].data

    # extract basic info
    vizier_epochs = EPOCH_CONFIG._get_section_by_query_kind("vizier")
    gaia_epoch = vizier_epochs["gaia"]
    non_gaia_epoch = vizier_epochs[survey]
    lat_col = survey_info["lat_column"]
    lon_col = survey_info["lon_column"]

    # return uncorrected detections
    if disable_corrections:
        gaia_data = pd.DataFrame()

    non_gaia_coords = SkyCoord(
        ra=non_gaia_data[lon_col].to_numpy() * u.deg,
        dec=non_gaia_data[lat_col].to_numpy() * u.deg,
        frame=survey_info["frame"],
        obstime=non_gaia_epoch,
        pm_ra_cosdec=np.zeros(len(non_gaia_data)) * u.mas / u.yr,
        pm_dec=np.zeros(len(non_gaia_data)) * u.mas / u.yr,
    )

    if not gaia_data.empty:
        # keep only necessary gaia columns
        gaia_data = gaia_data[["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"]]

        # discard any gaia sources that have nan positions/proper motions
        gaia_data = gaia_data.dropna()

        # combine non_gaia Vizier parameters into a single array
        params = []
        for key, val in survey_info.items():
            if key == "frame":
                continue
            if not isinstance(val, list):
                val = [val]
            params += val

        # keep only necessary columns
        non_gaia_data = non_gaia_data[params]

        # set up gaia and non-gaia arrays of SkyCoords
        gaia_coords = SkyCoord(
            ra=gaia_data["RA_ICRS"].to_numpy() * u.deg,
            dec=gaia_data["DE_ICRS"].to_numpy() * u.deg,
            pm_ra_cosdec=gaia_data["pmRA"].to_numpy() * (u.mas / u.yr),
            pm_dec=gaia_data["pmDE"].to_numpy() * (u.mas / u.yr),
            frame="icrs",
            obstime=gaia_epoch,
        )

        # match frames and correct gaia detections to non-gaia epoch
        non_gaia_coords = non_gaia_coords.transform_to(gaia_coords.frame)
        gaia_coords = gaia_coords.apply_space_motion(non_gaia_epoch)

        # nearest-neighbour cross match + check for those within radius of gaia sources
        index, separation, _ = match_coordinates_sky(non_gaia_coords, gaia_coords)

        # mask contains rows
        mask = separation < piggyback_radius * u.arcsec
        matched_index = np.where(mask)[0]

        # Prepare arrays
        pm_ra = np.zeros(len(non_gaia_coords)) * u.mas / u.yr
        pm_dec = np.zeros(len(non_gaia_coords)) * u.mas / u.yr

        # Fill matched entries
        pm_ra[matched_index] = gaia_coords.pm_ra_cosdec[index[matched_index]]
        pm_dec[matched_index] = gaia_coords.pm_dec[index[matched_index]]

        # Create a new SkyCoord in Gaia frame
        non_gaia_coords = SkyCoord(
            ra=non_gaia_coords.ra,
            dec=non_gaia_coords.dec,
            frame=gaia_coords.frame,
            obstime=non_gaia_coords.obstime,
            pm_ra_cosdec=pm_ra,
            pm_dec=pm_dec,
        )

        # correct all detections with proper motion information to the epoch of the image
        non_gaia_coords = non_gaia_coords.apply_space_motion(image.epoch)
    else:
        mask = [False] * len(non_gaia_data)

    # set up overlay DataFrame
    df = pd.DataFrame()
    df["survey"] = [survey] * len(non_gaia_coords)
    df["ra"] = non_gaia_coords.ra.deg
    df["dec"] = non_gaia_coords.dec.deg
    df["pm_ra_cosdec"] = non_gaia_coords.pm_ra_cosdec.to(u.mas / u.yr).value
    df["pm_dec"] = non_gaia_coords.pm_dec.to(u.mas / u.yr).value
    df["gaia_match"] = mask

    # duplicate above DataFrame for each requested magnitude + fill in these columns
    per_mag_dfs = []
    for mag, err in zip(survey_info["mags"], survey_info["errors"]):
        df["mag_name"] = mag
        df["mag"] = non_gaia_data[mag]
        df["err_name"] = err
        df["err"] = non_gaia_data[err]

        per_mag_dfs.append(df.copy())
    final_df = pd.concat(per_mag_dfs).reset_index(drop=True)

    # replace zero proper motion back to nan
    final_df[["pm_ra_cosdec", "pm_dec"]] = final_df[["pm_ra_cosdec", "pm_dec"]].replace(0, np.nan)

    # cull detections that are outside the final image bounds
    n_pixels = (image.hdu.data.shape[1], image.hdu.data.shape[0])
    pixel_scales = proj_plane_pixel_scales(image.wcs)
    x_bounds = (
        image.search_pos.ra.value - n_pixels[0] / 2 * pixel_scales[0],
        image.search_pos.ra.value + n_pixels[0] / 2 * pixel_scales[0],
    )
    y_bounds = (
        image.search_pos.dec.value - n_pixels[1] / 2 * pixel_scales[1],
        image.search_pos.dec.value + n_pixels[1] / 2 * pixel_scales[1],
    )
    ra_mask = (final_df["ra"] < x_bounds[0]) | (final_df["ra"] > x_bounds[1])
    dec_mask = (final_df["dec"] < y_bounds[0]) | (final_df["dec"] > y_bounds[1])
    cull_mask = ra_mask | dec_mask
    final_df = final_df.drop(final_df[cull_mask].index)

    if final_df.empty:
        return final_df

    # correct all detections with proper motion information to J2000 and search for SIMBAD IDs
    coord = dataframe_to_skycoord(final_df, image.epoch)
    coord = correct_skycoord(coord, image.epoch, Time("2000-01-01", format="iso"))

    # get SIMBAD object IDs
    ids = get_ids(coord, BASE_CONFIG._get("overlay_settings", "simbad_radius") * u.arcsec)

    if ids is RETURNS.EXCEPTION:
        return ids
    elif ids is RETURNS.NULL:
        final_df["simbad_id"] = "None"
    else:
        final_df["simbad_id"] = ids

    return final_df


def get_overlay(target: Target, image: Image, **kwargs: dict):
    """
    Fetches detection overlay information within a given image for a list of Vizier catalogue aliases or a dict of survey:band keys
    """

    overlay_dict = OVERLAY_CONFIG._as_dict()
    disable_corrections = kwargs.get("disable_corrections", False)

    overlays = kwargs.get("overlays")
    if not overlays:
        return None

    if isinstance(overlays, str):
        overlays = [overlays]

    # extract overlay info for requested surveys
    overlay_info = {}
    for survey in overlays:
        for section in overlay_dict.values():
            if survey in section:
                overlay_info[survey] = section[survey]
    if not overlay_info:
        raise ValueError("Failed to find any of the requested overlays in Overlay Config.")

    if isinstance(overlays, list):
        # keep only the first magnitude + error columns
        reduced_overlay_info = {}
        for survey, entry in overlay_info.items():
            new_entry = {}
            for key, val in entry.items():
                if isinstance(val, list):
                    new_entry[key] = val[:1]
                else:
                    new_entry[key] = val
            reduced_overlay_info[survey] = new_entry
        overlay_info = reduced_overlay_info

    elif isinstance(overlays, dict):
        reduced_overlay_info = {}
        for survey, entry in overlay_info.items():
            requested_mags, errors = [], []
            for mag in overlays[survey]:
                if mag not in entry["mags"]:
                    raise ValueError(f"Magnitude column '{mag}' not found in overlay definition for survey '{survey}'.")
                requested_mags.append(mag)
                errors.append(entry["errors"][entry["mags"].index(mag)])

            overlay_info[survey]["mags"] = requested_mags
            overlay_info[survey]["errors"] = errors

    overlay_data = []
    for survey, info in overlay_info.items():
        data = get_overlay_data(image, target, survey, info, disable_corrections)
        # if an exception is encountered, return EXCEPTION and set overlay=None, exception=True in structure
        if data is RETURNS.EXCEPTION:
            return data
        overlay_data.append(data)

    # combine overlay data from all requested surveys
    final_overlay = pd.concat(overlay_data).reset_index(drop=True)
    if final_overlay.empty:
        return None

    print(final_overlay)

    return final_overlay
