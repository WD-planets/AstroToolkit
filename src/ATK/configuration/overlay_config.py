from .alias_config import ALIAS_CONFIG
from .defaults.overlay_defaults import DEFAULTS, PATH
from .epoch_config import EPOCH_CONFIG
from .yaml_config.YAMLConfig import YAMLConfig

BASE_INDENT = 4

ALL_ALIASES = ALIAS_CONFIG._as_flattened_dict()
EPOCHS = EPOCH_CONFIG._as_dict()["vizier_aliases"]


def pprint_overlays(data, indent=0):
    """
    Pretty-print nested dicts of overlay data
    """

    lines = []

    for section, section_data in data.items():
        # Section header
        lines.append(f"[{section}]")

        for survey, survey_data in section_data.items():
            if survey in ALL_ALIASES:
                alias_str = ALL_ALIASES[survey]
            else:
                alias_str = "n/a"

            if survey in EPOCHS:
                epoch_str = EPOCHS[survey]
            else:
                epoch_str = "n/a"

            # Survey key: zero-indented
            lines.append(f"{' ' * int(BASE_INDENT / 2)}<{survey} ({alias_str} | {epoch_str})>")

            # Nested keys inside survey: indented
            for key, value in survey_data.items():
                indent = " " * BASE_INDENT
                if isinstance(value, list):
                    lines.append(f"{indent}{key}: {', '.join(map(str, value))}")
                else:
                    lines.append(f"{indent}{key} = {value}")

            # Blank line between surveys
            lines.append("")

    return "\n".join(lines)


class OverlayConfig(YAMLConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, None, None, pprint_overlays)


OVERLAY_CONFIG = OverlayConfig()
