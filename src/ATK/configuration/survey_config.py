import warnings

from astropy.time import Time

from .defaults.survey_defaults import DEFAULTS, PATH
from .yaml_config.YAMLConfig import YAMLConfig

BASE_INDENT = 4


def epoch_translator(value):
    """
    Translates strings to astropy ISOT Time objects
    """

    try:
        return Time(value, format="isot")
    except ValueError:
        warnings.warn("ATK: Invalid epoch found in config file. Use 'ATKepoch show' to see the invalid entry.")

        return "<Invalid ISOT Time Format>"


def pprint(data, indent=0):
    """
    Pretty-print nested dicts of overlay data
    """

    lines = []

    for section, section_data in data.items():
        # Section header
        lines.append(f"[{section}]")

        for survey, survey_data in section_data.items():
            # Survey key: zero-indented
            lines.append(f"{' ' * int(BASE_INDENT / 2)}<{survey}>")

            # Nested keys inside survey: indented
            for key, value in survey_data.items():
                indent = " " * BASE_INDENT
                if isinstance(value, list):
                    lines.append(f"{indent}{key}: {', '.join(map(str, value))}")
                else:
                    lines.append(f"{indent}{key} = {value}")

            # Blank line between surveys
            lines.append("")
        lines.append("")

    return "\n".join(lines)


class SurveyConfig(YAMLConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, None, None, pprint)

    def _get_aliases(self):
        self._load()

        aliases = {}
        for survey, survey_data in self._as_dict()["vizier"].items():
            if "id" in survey_data:
                aliases[survey] = survey_data["id"]

        return aliases

    def _get_epochs(self, kind: str):
        self._load()

        epochs = {}
        for survey, survey_data in self._as_dict()[kind].items():
            if "epoch" in survey_data:
                epochs[survey] = epoch_translator(survey_data["epoch"])

        return epochs

    def _get_overlays(self):
        self._load()

        positional_keys = ["lon_column", "lat_column", "frame"]
        photometric_keys = positional_keys + ["mags", "errors"]

        overlays = {"photometric": {}, "positional": {}}
        for survey, survey_data in self._as_dict()["vizier"].items():
            if all(key in survey_data.keys() for key in photometric_keys):
                overlays["photometric"][survey] = {key: survey_data[key] for key in photometric_keys}
            elif all(key in survey_data.keys() for key in positional_keys):
                overlays["positional"][survey] = {key: survey_data[key] for key in positional_keys}

        return overlays


SURVEY_CONFIG = SurveyConfig()
