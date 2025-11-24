import warnings

from astropy.time import Time

from .defaults.epoch_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig


def translator(value):
    try:
        return Time(value, format="iso")
    except ValueError:
        warnings.warn("ATK: Invalid epoch found in config file. Use 'ATKepoch show' to see the invalid entry.")
        return "<Invalid ISO Time Format>"


class EpochConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)


EPOCH_CONFIG = EpochConfig()
