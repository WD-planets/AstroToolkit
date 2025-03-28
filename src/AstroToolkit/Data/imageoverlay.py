import math

import numpy as np

from ..Configuration.baseconfig import ConfigStruct
from ..Configuration.epochs import EpochStruct
from ..Configuration.overlays import OverlayStruct
from ..Data.simbad_query import pos_query
from ..Misc.pmcorrection import correctradius
from ..PackageInfo import SurveyInfo
from ..Tools import correctpm, query

overlayInfo = OverlayStruct().read_overlays()
dataSurveyInfo = SurveyInfo().dataSurveyInfo

epochs = EpochStruct().epoch_list

config = ConfigStruct()
config.read_config()

piggyback_radius = int(config.overlay_piggyback_radius)


class OverlayData(object):
    def __init__(self, returned_data, survey):
        self.returned_data = returned_data
        self.survey = survey
        self.corrected_systems = []

    def correct_gaia_to_non_gaia(self):
        for i in range(0, len(self.returned_data["gaia"]["sources"])):
            if not math.isnan(self.returned_data["gaia"]["pmra"][i]) and not math.isnan(
                self.returned_data["gaia"]["pmdec"][i]
            ):
                (
                    self.returned_data["gaia"]["ra"][i],
                    self.returned_data["gaia"]["dec"][i],
                ) = correctpm(
                    input_time=epochs["gaia"],
                    target_time=epochs[self.survey],
                    pos=[
                        self.returned_data["gaia"]["ra"][i],
                        self.returned_data["gaia"]["dec"][i],
                    ],
                    pmra=self.returned_data["gaia"]["pmra"][i],
                    pmdec=self.returned_data["gaia"]["pmdec"][i],
                )

    def do_piggyback_correction(self):
        for i in range(0, len(self.returned_data["non_gaia"]["ra"])):
            for j in range(0, len(self.returned_data["gaia"]["sources"])):
                delta = (
                    math.sqrt(
                        (
                            self.returned_data["non_gaia"]["ra"][i]
                            - self.returned_data["gaia"]["ra"][j]
                        )
                        ** 2
                        + (
                            self.returned_data["non_gaia"]["dec"][i]
                            - self.returned_data["gaia"]["dec"][j]
                        )
                        ** 2
                    )
                    * 3600
                )
                if delta < piggyback_radius:
                    if not math.isnan(
                        self.returned_data["gaia"]["pmra"][j]
                    ) and not math.isnan(self.returned_data["gaia"]["pmdec"][j]):
                        (
                            self.returned_data["non_gaia"]["ra"][i],
                            self.returned_data["non_gaia"]["dec"][i],
                        ) = correctpm(
                            input_time=epochs[self.survey],
                            target_time=image_time,
                            pos=[
                                self.returned_data["non_gaia"]["ra"][i],
                                self.returned_data["non_gaia"]["dec"][i],
                            ],
                            pmra=self.returned_data["gaia"]["pmra"][j],
                            pmdec=self.returned_data["gaia"]["pmdec"][j],
                        )
                        self.returned_data["non_gaia"]["correction_pmra"][i] = (
                            self.returned_data["gaia"]["pmra"][j]
                        )
                        self.returned_data["non_gaia"]["correction_pmdec"][i] = (
                            self.returned_data["gaia"]["pmdec"][j]
                        )
                        self.corrected_systems.append(i)

    def scale_magnitudes(self):
        for mag in overlayInfo[self.survey]["mag_names"]:
            self.returned_data["non_gaia"][f"{mag}_marker_size"] = []
            for i in range(0, len(self.returned_data["non_gaia"]["ra"])):
                if not math.isnan(self.returned_data["non_gaia"][mag][i]):
                    radius_multiplier = self.returned_data["non_gaia"][mag][i] / 20.7
                    base_marker_size = 0.75
                    self.returned_data["non_gaia"][f"{mag}_marker_size"].append(
                        (
                            half_image_size / 50
                            + (half_image_size / 75) ** (radius_multiplier * 1.15)
                        )
                        * base_marker_size
                        + 0.0005
                    )
                else:
                    self.returned_data["non_gaia"][f"{mag}_marker_size"].append(None)


def get_overlay_data(data, survey):
    params = overlayInfo[survey]
    if params["overlay_type"] != "tracer":
        obj_id_name = overlayInfo[survey]["id_name"]

    if data.source:
        radius = correctradius(
            source=data.source,
            input_time=epochs["gaia"],
            target_time=epochs[survey],
            radius=data.data["size"],
        )
    else:
        radius = data.data["size"]

    if (
        params["overlay_type"] == "detection_mag"
        or params["overlay_type"] == "detection"
    ):
        non_gaia_systems = query(
            kind="data",
            survey=survey,
            pos=data.data["image_focus"],
            radius=radius,
            level="internal",
        ).data
        if not non_gaia_systems:
            return None
        gaia_systems = query(
            kind="data",
            survey="gaia",
            pos=data.data["image_focus"],
            radius=radius,
            level="internal",
        ).data

        # When generating a Gaia overlay, Gaia is considered a non_gaia survey. Correcting this using Gaia obviously does nothing
        returned_data = {
            "gaia": {
                "sources": gaia_systems["source_id"],
                "ra": gaia_systems["ra"],
                "dec": gaia_systems["dec"],
                "pmra": gaia_systems["pmra"],
                "pmdec": gaia_systems["pmdec"],
            },
            "non_gaia": {
                "ra": non_gaia_systems[params["ra_name"]],
                "dec": non_gaia_systems[params["dec_name"]],
                "correction_pmra": [
                    np.nan for x in non_gaia_systems[params["ra_name"]]
                ],
                "correction_pmdec": [
                    np.nan for x in non_gaia_systems[params["ra_name"]]
                ],
                "obj_id": non_gaia_systems[obj_id_name],
            },
        }

    if params["overlay_type"] == "detection_mag":
        for mag in overlayInfo["gaia"]["mag_names"]:
            returned_data["gaia"][mag] = gaia_systems[mag]

        for mag in params["mag_names"]:
            returned_data["non_gaia"][mag] = non_gaia_systems[mag]

    if (
        params["overlay_type"] == "detection_mag"
        or params["overlay_type"] == "detection"
    ):
        global image_time
        image_time = data.data["image_time"]
        global half_image_size
        half_image_size = data.data["size"] / 7200

        overlay_data = OverlayData(returned_data, survey)
        overlay_data.correct_gaia_to_non_gaia()
        overlay_data.do_piggyback_correction()
        if params["overlay_type"] == "detection_mag":
            overlay_data.scale_magnitudes()
        formatted_data, corrected_systems = (
            overlay_data.returned_data,
            overlay_data.corrected_systems,
        )

    if params["overlay_type"] == "detection_mag":
        overlay = []
        for mag in params["mag_names"]:
            for i in range(0, len(formatted_data["non_gaia"]["ra"])):
                if not math.isnan(formatted_data["non_gaia"][mag][i]):
                    overlay_entry = {
                        "overlay_type": params["overlay_type"],
                        "marker_type": params["marker_type"],
                        "corrected": i in corrected_systems,
                        "ra": formatted_data["non_gaia"]["ra"][i],
                        "dec": formatted_data["non_gaia"]["dec"][i],
                        "marker_size": formatted_data["non_gaia"][f"{mag}_marker_size"][
                            i
                        ],
                        "colour_index": params["colour_index"][
                            params["mag_names"].index(mag)
                        ],
                        "mag_name": mag,
                        "mag": formatted_data["non_gaia"][mag][i],
                        "survey": survey,
                        "obj_id": str(formatted_data["non_gaia"]["obj_id"][i]),
                        "correction_pmra": formatted_data["non_gaia"][
                            "correction_pmra"
                        ][i],
                        "correction_pmdec": formatted_data["non_gaia"][
                            "correction_pmdec"
                        ][i],
                    }
                    overlay.append(overlay_entry)
    elif params["overlay_type"] == "detection":
        overlay = []
        for i in range(0, len(formatted_data["non_gaia"]["ra"])):
            overlay_entry = {
                "survey": survey,
                "obj_id": str(formatted_data["non_gaia"]["obj_id"][i]),
                "overlay_type": params["overlay_type"],
                "marker_type": params["marker_type"],
                "marker_size": "NA",
                "corrected": i in corrected_systems,
                "ra": formatted_data["non_gaia"]["ra"][i],
                "dec": formatted_data["non_gaia"]["dec"][i],
                "j2000_ra": "",
                "j2000_dec": "",
                "colour_index": params["colour_index"],
                "mag_name": "NA",
                "correction_pmra": formatted_data["non_gaia"]["correction_pmra"][i],
                "correction_pmdec": formatted_data["non_gaia"]["correction_pmdec"][i],
            }
            overlay.append(overlay_entry)

    if params["overlay_type"] == "tracer":
        lightcurve_data = query(
            kind="lightcurve",
            survey=survey,
            pos=data.data["image_focus"],
            radius=radius,
            level="internal",
        ).data

        if not lightcurve_data:
            return None

        data_exists = False
        for band in lightcurve_data:
            if band["mag"]:
                data_exists = True
        if not data_exists:
            return None

        combined_ra, combined_dec = [], []
        for band in lightcurve_data:
            if band["ra"] and band["dec"]:
                combined_ra += band["ra"]
                combined_dec += band["dec"]

        overlay = []
        for ra, dec in zip(combined_ra, combined_dec):
            overlay_entry = {
                "survey": survey,
                "overlay_type": params["overlay_type"],
                "marker_type": params["marker_type"],
                "marker_size": "NA",
                "ra": ra,
                "dec": dec,
                "colour": params["colour_index"],
            }
            overlay.append(overlay_entry)

    return overlay


def overlay_query(data, overlays):
    overlays_data = []

    from ..Configuration.baseconfig import ConfigStruct

    config = ConfigStruct()
    config.read_config()

    if isinstance(overlays, list):
        for survey in overlays:
            if survey:
                if overlayInfo[survey]["overlay_type"] == "detection_mag":
                    overlayInfo[survey]["mag_names"] = [
                        overlayInfo[survey]["default_mag"]
                    ]
                overlay_data = get_overlay_data(data, survey)
                if overlay_data:
                    overlays_data += overlay_data
    elif isinstance(overlays, dict):
        for survey, mag_names in overlays.items():
            if mag_names:
                if not isinstance(mag_names, list):
                    mag_names = [mag_names]
                overlayInfo[survey]["mag_names"] = mag_names
            overlay_data = get_overlay_data(data, survey)
            if overlay_data:
                overlays_data += overlay_data

    if overlays_data == []:
        overlays_data = None

    for index, data_point in enumerate(overlays_data):
        if data_point["overlay_type"] != "tracer":
            if data_point["corrected"]:
                search_ra, search_dec = correctpm(
                    input_time=image_time,
                    target_time=[2000, 0],
                    pos=[data_point["ra"], data_point["dec"]],
                    pmra=data_point["correction_pmra"],
                    pmdec=data_point["correction_pmdec"],
                )
            else:
                search_ra, search_dec = data_point["ra"], data_point["dec"]
            identifier = pos_query([search_ra, search_dec])
            data_point["simbad_id"] = str(identifier)
        overlays_data[index] = data_point

    # cull data points that are outside the image (also need to actually fix the query radius at some point)
    image_header = data.data["image_header"]
    xlim, ylim = image_header["NAXIS1"], image_header["NAXIS2"]

    if xlim != ylim:
        xlim = ylim

    wcs = data.data["wcs"]
    x_points, y_points = (
        np.arange(start=0, stop=xlim + 1, step=1),
        np.arange(start=0, stop=ylim + 1, step=1),
    )

    coords = wcs.all_pix2world(x_points, y_points, 1)
    x_points, y_points = coords[0], coords[1]

    bad_indices = []
    for index, data_point in enumerate(overlays_data):
        if not (min(x_points) <= data_point["ra"] <= max(x_points)):
            bad_indices.append(index)
            continue
        if not (min(y_points) <= data_point["dec"] <= max(y_points)):
            bad_indices.append(index)
            continue

    overlays_data = [val for i, val in enumerate(overlays_data) if i not in bad_indices]

    return overlays_data
