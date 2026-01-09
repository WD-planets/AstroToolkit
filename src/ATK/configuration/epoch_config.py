import warnings

from astropy.time import Time
from erfa import ErfaWarning

from .defaults.epoch_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig

# ignore bad distance warning
warnings.filterwarnings("ignore", category=ErfaWarning)


def translator(value):
    """
    Translates strings to astropy ISOT Time objects
    """

    try:
        return Time(value, format="isot")
    except ValueError:
        warnings.warn("ATK: Invalid epoch found in config file. Use 'ATKepoch show' to see the invalid entry.")

        return "<Invalid ISOT Time Format>"


class EpochConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)

    def get_section_by_query_kind(self, query_kind: str):
        """
        Returns a requested epoch config section (i.e. the section for a specified query kind) as a dict
        """

        if query_kind == "vizier":
            section = f"{query_kind}_aliases"
        else:
            section = f"{query_kind}_surveys"

        return self.as_dict()[section]


EPOCH_CONFIG = EpochConfig()
