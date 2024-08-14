import math

from ..Data.dataquery import SurveyInfo
from ..Misc.pmcorrection import correctradius
from ..Tools import correctpm, query

survey_times = SurveyInfo().times
survey_params = SurveyInfo().overlay_param_names

from ..Configuration.baseconfig import ConfigStruct

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
                    input_time=survey_times["gaia"],
                    target_time=survey_times[self.survey],
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
                            input_time=survey_times[self.survey],
                            target_time=image_time,
                            pos=[
                                self.returned_data["non_gaia"]["ra"][i],
                                self.returned_data["non_gaia"]["dec"][i],
                            ],
                            pmra=self.returned_data["gaia"]["pmra"][j],
                            pmdec=self.returned_data["gaia"]["pmdec"][j],
                        )
                        self.corrected_systems.append(i)

    def scale_magnitudes(self):
        for mag in survey_params[self.survey]["mag_names"]:
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
    params = survey_params[survey]

    if data.source:
        radius = correctradius(
            source=data.source,
            input_time=survey_times["gaia"],
            target_time=survey_times[survey],
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
            },
        }

    if params["overlay_type"] == "detection_mag":
        for mag in survey_params["gaia"]["mag_names"]:
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
                    }
                    overlay_entry["marker_size"] = formatted_data["non_gaia"][
                        f"{mag}_marker_size"
                    ][i]
                    overlay_entry["colour"] = params["colours"][
                        params["mag_names"].index(mag)
                    ]
                    overlay_entry["mag_name"] = mag
                    overlay_entry["survey"] = survey
                    overlay.append(overlay_entry)
    elif params["overlay_type"] == "detection":
        overlay = []
        for i in range(0, len(formatted_data["non_gaia"]["ra"])):
            overlay_entry = {
                "survey": survey,
                "overlay_type": params["overlay_type"],
                "marker_type": params["marker_type"],
                "corrected": i in corrected_systems,
                "ra": formatted_data["non_gaia"]["ra"][i],
                "dec": formatted_data["non_gaia"]["dec"][i],
                "colour": params["colour"],
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
                "ra": ra,
                "dec": dec,
                "colour": params["colour"],
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
                if survey_params[survey]["overlay_type"] == "detection_mag":
                    survey_params[survey]["mag_names"] = [
                        getattr(config, f"{survey}_overlay_mag")
                    ]
                overlay_data = get_overlay_data(data, survey)
                if overlay_data:
                    overlays_data += overlay_data
    elif isinstance(overlays, dict):
        for survey, mag_names in overlays.items():
            if mag_names:
                if not isinstance(mag_names, list):
                    mag_names = [mag_names]
                survey_params[survey]["mag_names"] = mag_names
            overlay_data = get_overlay_data(data, survey)
            if overlay_data:
                overlays_data += overlay_data

    if overlays_data == []:
        overlays_data = None

    return overlays_data
