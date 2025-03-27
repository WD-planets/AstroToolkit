import os

import yaml

from ..PackageInfo import SurveyInfo

packageInfo = SurveyInfo()
overlayInfo = packageInfo.overlay_param_names
magInfo = {}
for survey in packageInfo.sed_param_names:
    magInfo[survey] = packageInfo.sed_param_names[survey]["mag_names"]


class overlayStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.default_surveys = packageInfo.supported_overlays

        self.overlay_file = files("AstroToolkit.Configuration").joinpath("ATKoverlays.ini")
        if not os.path.isfile(self.overlay_file):
            print("No ATKoverlays.ini found. Generating one with default values...")
            self.default_setup()

    def default_setup(self):
        defaults = {"scaled_detections", "detections", "tracers"}
        for survey in self.default_surveys:
            defaults[survey] = {}
            mag_str = ""
            for mag in magInfo[survey]:
                mag_str += mag
            defaults[survey]["magnitudes"] = mag_str

        with open(self.overlay_file, "w") as file:
            yaml.dump(defaults, file)
