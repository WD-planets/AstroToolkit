from configparser import ConfigParser
from pathlib import Path
from types import FunctionType

from ...utilities.file_io import open_file


class ConfigSection:
    """
    Proxy object for a config section
    """

    def __init__(self, parent, name):
        self._parent = parent
        self._name = name

    def __getattr__(self, key):
        config = self._parent._(self._name, {})
        if key not in config:
            raise AttributeError(f"Key '{key}' not found in section '{self._name}'")
        return config[key]

    def __setattr__(self, key, value):
        if key in {"_parent", "_name"}:
            super().__setattr__(key, value)
            return
        self._parent._set(self._name, key, value)

    def __delattr__(self, key):
        self._parent._del(self._name, key)


def parser_to_dict(parser: ConfigParser, translator: FunctionType) -> dict:
    """
    Convert a ConfigParser to a dictionary of values, using a translation function if provided.
    """

    config = {}
    for section in parser.sections():
        config[section] = {key: translator(val) for key, val in parser[section].items()}

    return config


class ParserConfig:
    """
    Base class for a ConfigParser/dict-style configuration manager
    """

    def __init__(self, path: Path, defaults: dict, translator: FunctionType | None = None):
        self._path = path
        self._defaults = defaults
        self._translator = translator or (lambda x: x)

        self._parser = None
        self._config = None
        self._load()

    def _load(self) -> None:
        """
        Load parser and config from file if it exists, if not create one with default values
        """

        parser = ConfigParser()

        if self._path.exists():
            parser.read(self._path)
        else:
            parser.read_dict(self._defaults)
            self._parser = parser
            self._save()

        self._config = parser_to_dict(parser, self._translator)
        self._parser = parser

    def _save(self) -> None:
        """
        Save parser to file
        """

        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "w") as f:
            self._parser.write(f)

    def reset(self) -> None:
        """
        Reset parser and config to their default states, and write this to the config file
        """

        parser = ConfigParser()
        parser.read_dict(self._defaults)

        with open(self._path, "w") as f:
            parser.write(f)
        self._load()

    def _set(self, section: str, key: str, value: any, only_existing=False) -> None:
        """
        Set the value of a given config key in a given section
        """

        self._load()
        if section not in self._parser:
            raise ValueError(f"Section '{section}' not found in config.")
        if only_existing and key not in self._parser[section]:
            raise ValueError(f"Key '{key}' not found in section '{section}' of config.")
        self._parser.set(section, key, str(value))
        self._config[section][key] = value
        self._save()

    def _del(self, section: str, key: str) -> None:
        """
        Deletes a given key in a given section of the config
        """
        self._load()
        if section not in self._parser:
            raise ValueError(f"Section '{section}' not found in config.")
        if key not in self._parser[section]:
            raise ValueError(f"Key '{key}' not found in section '{section}'.")
        del self._config[section][key]
        del self._parser[section][key]
        self._save()

    def _get(self, section: str, key: str, fallback=None):
        """
        Get the value of a key in a given section
        """

        self._load()
        return self._config.get(section, {}).get(key, fallback)

    def _get_section(self, section: str):
        """
        Get an entire section by its key in the config
        """

        self._load()
        return self._config.get(section, {})

    def _as_dict(self):
        """
        Get entire config as a 2D dictionary
        """

        self._load()
        return self._config

    def show(self) -> None:
        """
        Print config to stdout
        """

        self._load()
        for section, key_vals in self._config.items():
            print(f"[{section}]")
            for key, val in key_vals.items():
                print(f"{key} = {val}")
            print()

    def open(self) -> None:
        """
        Opens the config file in the default text editor
        """

        open_file(self._path)

    def __getattr__(self, section):
        if section in self.__dict__:
            return self.__dict__[section]

        if section not in self._config:
            raise ValueError(f"Config section '{section}' does not exist.")

        return ConfigSection(self, section)
