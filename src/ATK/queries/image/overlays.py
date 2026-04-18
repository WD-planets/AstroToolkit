import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs.utils import proj_plane_pixel_scales

from ...configuration.base_config import BASE_CONFIG
from ...configuration.survey_config import SURVEY_CONFIG
from ...structures.Image import Image
from ...structures.Target import Target
from ...Tools.query import query
from ...utilities.coordinates import correct_coords, correct_radius
from ...utilities.defaults import RETURNS
from ..simbad.simbad_query import get_ids


def check_finite(arr):
    return np.isfinite(arr) & (arr is not None)


def _get_unit(column):
    return getattr(column, "unit", None) or u.deg


def get_overlay_data(image: Image, target: int | SkyCoord, survey: str, survey_info: dict, disable_corrections=False) -> pd.DataFrame:
    radius = correct_radius(target, image.size, "vizier", survey) * 1.25
    piggyback_radius = BASE_CONFIG._get("overlay_settings", "piggyback_radius")

    non_gaia_data = query(kind="vizier", targets=image.search_pos, radius=radius, survey=survey).data

    if not non_gaia_data:
        return pd.DataFrame()
    non_gaia_data = non_gaia_data[0].table

    gaia_data = query(kind="vizier", targets=image.search_pos, radius=radius, survey="gaia").data

    if not gaia_data:
        gaia_data = pd.DataFrame()
    else:
        gaia_data = gaia_data[0].table.to_pandas()

    vizier_epochs = SURVEY_CONFIG._get_epochs("vizier")
    gaia_epoch = vizier_epochs["gaia"]
    non_gaia_epoch = vizier_epochs[survey]

    lat_col = survey_info["lat"]
    lon_col = survey_info["lon"]

    ra_unit = _get_unit(non_gaia_data[lon_col])
    dec_unit = _get_unit(non_gaia_data[lat_col])

    if disable_corrections:
        gaia_data = Table()

    non_gaia_coords = [
        SkyCoord(
            ra=row[lon_col] * ra_unit,
            dec=row[lat_col] * dec_unit,
            frame=survey_info["frame"],
            pm_ra_cosdec=np.nan * u.mas / u.yr,
            pm_dec=np.nan * u.mas / u.yr,
            obstime=non_gaia_epoch,
        )
        for row in non_gaia_data
    ]

    if len(gaia_data):
        gaia_data = gaia_data[["RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Plx"]]

        mask = check_finite(gaia_data["Plx"]) & (gaia_data["Plx"] > 0)
        gaia_data["dist"] = np.where(mask, 1000 / gaia_data["Plx"], np.nan)

        gaia_data = gaia_data.dropna(subset=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE"])

        params = []
        for key, val in survey_info.items():
            if key == "frame":
                continue
            if not isinstance(val, list):
                val = [val]
            params += val

        non_gaia_data = non_gaia_data[params]

        gaia_coords = []
        for row in gaia_data.itertuples(index=False):
            ra = row.RA_ICRS * u.deg
            dec = row.DE_ICRS * u.deg
            pmra = row.pmRA * (u.mas / u.yr)
            pmdec = row.pmDE * (u.mas / u.yr)

            dist_val = row.dist * u.pc if np.isfinite(row.dist) else None
            if isinstance(dist_val, Quantity) or dist_val is None:
                pass
            elif np.isnan(dist_val):
                raise ValueError("Bad distance.")

            # test for bad distance propagation
            # rng = np.random.uniform(0, 1)
            # if rng < 0.5:
            #     dist_val = None

            gaia_coords.append(SkyCoord(ra=ra, dec=dec, frame="icrs", pm_ra_cosdec=pmra, pm_dec=pmdec, distance=dist_val, obstime=gaia_epoch))

        non_gaia_coords = [c.transform_to("icrs") for c in non_gaia_coords]
        gaia_coords = correct_coords(gaia_coords, non_gaia_epoch)

        index = np.empty(len(non_gaia_coords), dtype=int)
        separation = np.empty(len(non_gaia_coords)) * u.arcsec

        gaia_stack = SkyCoord(ra=[c.ra for c in gaia_coords], dec=[c.dec for c in gaia_coords], frame="icrs")

        for i, src in enumerate(non_gaia_coords):
            sep = src.separation(gaia_stack)
            j = np.argmin(sep)

            index[i] = j
            separation[i] = sep[j]

        mask = separation < piggyback_radius * u.arcsec
        matched_index = np.where(mask)[0]

        pm_ra = np.full(len(non_gaia_coords), np.nan)
        pm_dec = np.full(len(non_gaia_coords), np.nan)
        dist = np.full(len(non_gaia_coords), np.nan)

        if len(matched_index):
            pm_ra[matched_index] = [gaia_coords[i].pm_ra_cosdec.to_value(u.mas / u.yr) for i in index[matched_index]]
            pm_dec[matched_index] = [gaia_coords[i].pm_dec.to_value(u.mas / u.yr) for i in index[matched_index]]
            dist[matched_index] = [gaia_coords[i].distance.to_value(u.pc) if gaia_coords[i].distance.unit is not u.one else np.nan for i in index[matched_index]]

        new_coords = []

        for i, c in enumerate(non_gaia_coords):
            dist_val = dist[i]
            pmra_val = pm_ra[i] * u.mas / u.yr
            pmdec_val = pm_dec[i] * u.mas / u.yr

            dist_val = dist[i] * u.pc if np.isfinite(dist[i]) else None
            if isinstance(dist_val, Quantity) or dist_val is None:
                pass
            elif np.isnan(dist_val):
                raise ValueError("Bad distance.")

            new_coords.append(SkyCoord(ra=c.ra, dec=c.dec, frame="icrs", obstime=c.obstime, pm_ra_cosdec=pmra_val, pm_dec=pmdec_val, distance=dist_val))

        non_gaia_coords = new_coords
        non_gaia_coords, correction = correct_coords(non_gaia_coords, image.epoch, get_correction=True)

    else:
        correction = ["none"] * len(non_gaia_data)

    df = pd.DataFrame(
        {
            "survey": [survey] * len(non_gaia_coords),
            "ra": [c.ra.deg for c in non_gaia_coords],
            "dec": [c.dec.deg for c in non_gaia_coords],
            "pm_ra_cosdec": [c.pm_ra_cosdec.to_value(u.mas / u.yr) for c in non_gaia_coords],
            "pm_dec": [c.pm_dec.to_value(u.mas / u.yr) for c in non_gaia_coords],
            "dist": [c.distance.to_value(u.pc) if c.distance.unit is not u.one else np.nan for c in non_gaia_coords],
            "gaia_match": correction,
        }
    )

    if "mags" in survey_info and "errors" in survey_info:
        per_mag_dfs = []

        for mag, err in zip(survey_info["mags"], survey_info["errors"]):
            tmp = df.copy()
            tmp["mag_name"] = mag
            tmp["mag"] = np.asarray(non_gaia_data[mag])
            tmp["err_name"] = err
            tmp["err"] = np.asarray(non_gaia_data[err])

            per_mag_dfs.append(tmp)

        final_df = pd.concat(per_mag_dfs).reset_index(drop=True)
    else:
        df["mag_name"] = "n/a"
        df["mag"] = np.nan
        df["err"] = "n/a"
        df["err_name"] = np.nan
        final_df = df

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

    final_df = final_df.drop(final_df[(ra_mask | dec_mask)].index)

    if final_df.empty:
        return final_df

    final_coords = []

    for _, row in final_df.iterrows():
        dist_val = row["dist"] * u.pc if np.isfinite(row["dist"]) else None

        if isinstance(dist_val, Quantity) or dist_val is None:
            pass
        elif np.isnan(dist_val):
            raise ValueError("Bad distance.")

        final_coords.append(
            SkyCoord(
                ra=row["ra"] * u.deg,
                dec=row["dec"] * u.deg,
                pm_ra_cosdec=row["pm_ra_cosdec"] * (u.mas / u.yr),
                pm_dec=row["pm_dec"] * (u.mas / u.yr),
                distance=dist_val,
                frame="icrs",
                obstime=image.epoch,
            )
        )

    final_coords = correct_coords(final_coords, Time("2000-01-01", format="iso"))

    ids = get_ids(final_coords, BASE_CONFIG._get("overlay_settings", "simbad_radius") * u.arcsec)

    if ids is RETURNS.EXCEPTION:
        return ids
    elif ids is RETURNS.NULL:
        final_df["simbad_id"] = "None"
    else:
        final_df["simbad_id"] = ids

    final_df = final_df.rename(columns={"gaia_match": "correction"})

    return final_df


def get_overlay(target: Target, image: Image, **kwargs: dict):
    """
    Fetches detection overlay information within a given image for a list of Vizier catalogue aliases or a dict of survey:band keys
    """

    overlay_dict = SURVEY_CONFIG._get_overlays()
    disable_corrections = kwargs.get("disable_correction", False)

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

    return final_overlay
