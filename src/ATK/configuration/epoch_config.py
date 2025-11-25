import warnings

from astropy.time import Time
from erfa import ErfaWarning

from .defaults.epoch_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig

# ignore bad distance warning
warnings.filterwarnings("ignore", category=ErfaWarning)


def translator(value):
    try:
        return Time(value, format="isot")
    except ValueError:
        print(value)
        warnings.warn("ATK: Invalid epoch found in config file. Use 'ATKepoch show' to see the invalid entry.")
        return "<Invalid ISOT Time Format>"


class EpochConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)

    def get_section_by_query_kind(self, query_kind: str):
        match query_kind:
            case "vizier":
                section = "vizier_aliases"
            case "lightcurve" | "spectrum":
                section = f"{query_kind}_surveys"

        return self.as_dict()[section]


EPOCH_CONFIG = EpochConfig()
