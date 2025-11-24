from astropy.time import Time

from .defaults.epoch_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig


def translator(value):
    return Time(value, format="iso")


class AliasConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)


EPOCH_CONFIG = AliasConfig()
