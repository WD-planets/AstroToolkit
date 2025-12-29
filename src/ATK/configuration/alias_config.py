from .defaults.alias_defaults import DEFAULTS, PATH
from .parser_config.ParserConfig import ParserConfig


class AliasConfig(ParserConfig):
    def __init__(self):
        super().__init__(PATH, DEFAULTS, None)

    def as_flattened_dict(self) -> None:
        """
        Returns config as a single flattened dictionary of alias:id definitions
        """

        data = {}
        for key, sub_dict in self._config.items():
            for key, val in sub_dict.items():
                data[key] = val

        return data


ALIAS_CONFIG = AliasConfig()
