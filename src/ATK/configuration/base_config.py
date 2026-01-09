from .defaults.base_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig


def translator(value):
    # bool
    if isinstance(value, bool):
        return value

    # string representation of None/bool
    if isinstance(value, str) and value.lower() == "true":
        return True
    elif isinstance(value, str) and value.lower() == "false":
        return False
    elif isinstance(value, str) and value.lower() == "none":
        return None

    # floats/integers
    try:
        val = float(value)
        if val.is_integer():
            return int(val)
        else:
            return float(value)
    except (ValueError, TypeError):
        return value


class BaseConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, translator)

    def _set(self, section, key, value):
        super()._set(section, key, value, True)


BASE_CONFIG = BaseConfig()
