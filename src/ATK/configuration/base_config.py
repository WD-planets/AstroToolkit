from .defaults.base_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig


def translator(value):
    match value:
        case "true":
            return True
        case "false":
            return False
        case "none":
            return None
        case _:
            pass

    try:
        return int(value)
    except ValueError:
        pass

    try:
        return float(value)
    except ValueError:
        return value


class BaseConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)


BASE_CONFIG = BaseConfig()
