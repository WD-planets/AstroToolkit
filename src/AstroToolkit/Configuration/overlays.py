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

        self.overlay_file = files("AstroToolkit.Configuration").joinpath("ATKoverlays.yaml")
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

        with open(self.overlay_file, "w") as file:
            yaml.dump(defaults, file, sort_keys=False, indent=4, Dumper=customDumper)

    def read_overlays(self):
        with open(self.overlay_file) as file:
            try:
                data = yaml.safe_load(file)
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
