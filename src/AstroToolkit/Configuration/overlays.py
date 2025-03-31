import os

import yaml

from ..PackageInfo import OverlayInfo, SurveyInfo

surveyInfo = SurveyInfo()
overlayInfo = OverlayInfo()
magInfo = {}


# overrides yaml dumper to print newlines between sections
class customDumper(yaml.SafeDumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()
            super().write_line_break()


class OverlayStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.overlay_file = files("AstroToolkit.Configuration").joinpath("ATKOverlays.yaml")
        if not os.path.isfile(self.overlay_file):
            print("No ATKoverlays.yaml found. Generating one with default values...")
            self.default_setup()

    def default_setup(self):
        defaults = {}

        scaled_detection_surveys = overlayInfo.detection_magSurveys
        for survey in scaled_detection_surveys:
            del scaled_detection_surveys[survey]["overlay_type"]

        detection_surveys = overlayInfo.detectionSurveys
        for survey in detection_surveys:
            del detection_surveys[survey]["overlay_type"]

        tracer_surveys = overlayInfo.tracerSurveys
        for survey in tracer_surveys:
            del tracer_surveys[survey]["overlay_type"]

        defaults["scaled_detections"] = scaled_detection_surveys
        defaults["detections"] = detection_surveys
        defaults["tracers"] = tracer_surveys

        # don't include tracers
        # del defaults["tracers"]

        with open(self.overlay_file, "w") as file:
            yaml.dump(defaults, file, sort_keys=False, indent=4, Dumper=customDumper)

    def read_overlays(self, raw=False):
        with open(self.overlay_file) as file:
            try:
                data = yaml.safe_load(file)
                if raw:
                    return data
            except yaml.YAMLError as e:
                print(e)

        overlay_dict = {}
        for survey, info in data["scaled_detections"].items():
            overlay_dict[survey] = info
            overlay_dict[survey]["overlay_type"] = "detection_mag"
            overlay_dict[survey]["marker_type"] = "circle"
            overlay_dict[survey]["default_mag"] = info["mag_names"][0]
        for survey, info in data["detections"].items():
            overlay_dict[survey] = info
            overlay_dict[survey]["overlay_type"] = "detection"
            overlay_dict[survey]["marker_type"] = "cross"
        for survey, info in data["tracers"].items():
            overlay_dict[survey] = info
            overlay_dict[survey]["overlay_type"] = "tracer"
            overlay_dict[survey]["marker_type"] = "cross"

        rolling_index = 0
        for survey, info in overlay_dict.items():
            indexes = []
            if info["overlay_type"] == "detection_mag":
                length = len(info["mag_names"])
            else:
                length = 1
            for i in range(0, length):
                indexes.append(rolling_index)
                rolling_index += 1
            overlay_dict[survey]["colour_index"] = indexes

        return overlay_dict

    @property
    def supportedOverlays(self):
        data = self.read_overlays()
        return list(data.keys())

    def add_overlay(self, survey, ra_name, dec_name, id_name, mag_names):
        from .catalogue_setup import CatalogueStruct

        aliases = CatalogueStruct().get_alias_list()
        if survey not in aliases:
            raise ValueError(
                f"No alias found for survey '{survey}'. One should first be added using: ATKalias add {survey} <Vizier ID>"
            )

        overlayData = self.read_overlays(raw=True)

        if mag_names:
            section = "scaled_detection"
            overlayData["scaled_detections"][survey] = {
                "ra_name": ra_name,
                "dec_name": dec_name,
                "id_name": id_name,
                "mag_names": mag_names,
            }
        else:
            section = "detection"
            overlayData["detections"][survey] = {"ra_name": ra_name, "dec_name": dec_name, "id_name": id_name}

        with open(self.overlay_file, "w") as file:
            yaml.dump(overlayData, file, sort_keys=False, indent=4, Dumper=customDumper)

        print(f"Added {section} overlay defintion for alias '{survey}'.")

    def del_overlay(self, section, survey):
        overlayData = self.read_overlays(raw=True)

        if survey in overlayData[f"{section}s"]:
            del overlayData[f"{section}s"][survey]
        else:
            raise ValueError(f"Could not find existing {section} overlay definition for alias '{survey}'.")

        with open(self.overlay_file, "w") as file:
            yaml.dump(overlayData, file, sort_keys=False, indent=4, Dumper=customDumper)

    def reset_overlays(self):
        self.default_setup()
