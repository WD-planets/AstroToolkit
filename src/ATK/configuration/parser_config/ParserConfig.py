from configparser import ConfigParser
from pathlib import Path
from types import FunctionType


def parser_to_dict(parser: ConfigParser, translator: FunctionType) -> dict:
    """
    Convert a ConfigParser to a dictionary of values, using a translation function if provided.
    """

    config = {}
    for section in parser:
        config[section] = {key: translator(val) for key, val in parser[section].items()}

    return config


class ParserConfig(object):
    """
    Base class for a ConfigParser/dict-style configuration manager
    """

    def __init__(self, path: Path, defaults: dict, translator: FunctionType | None = None):
        self._path = path
        self._defaults = defaults
        self._translator = translator or (lambda x: x)

        self._parser = None
        self._config = None
        self.load()

    def load(self) -> None:
        parser = ConfigParser()

        if self._path.exists():
            parser.read(self._path)
        else:
            parser.read_dict(self._defaults)

        self._config = parser_to_dict(parser, self._translator)
        self._parser = parser

    def save(self) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "w") as f:
            self._parser.write(f)

    def reset(self) -> None:
        parser = ConfigParser()
        parser.read_dict(self._defaults)

        with open(self._path) as f:
            parser.write(f)
        self.load()

    def set(self, section: str, key: str, value: any) -> None:
        if section not in self._parser:
            raise ValueError(f"Unexpected section '{section}'.")

        self._parser.set(section, key, str(value))
        self._config[section][key] = value

    def get(self, section: str, key: str, fallback=None):
        self.load()
        return self._config.get(section, {}).get(key, fallback)

    def get_section(self, section: str):
        return self._config.get(section, {})

    def as_dict(self):
        return self._config

    def show(self) -> None:
        for section, key_vals in self._config.items():
            print(f"[{section}]")
            for key, val in key_vals.items():
                print(f"{key} = {val}")
            print()
