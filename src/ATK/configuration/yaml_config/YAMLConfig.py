from pathlib import Path
from types import FunctionType

import yaml

from ...utilities.file_io import open_file
from .yaml_io import CustomDumper, default_printer

YAML_INDENT = 4


class _SectionProxy:
    """
    Acts as a live proxy to a section in a YAMLConfig
    """

    def __init__(self, parent: "YAMLConfig", section: str):
        self._parent = parent
        self._section = section

    def __getitem__(self, key):
        self._parent._load()
        if key not in self._parent._config[self._section]:
            self._parent._config[self._section][key] = None
        return _KeyProxy(self._parent, self._section, key)

    def __setitem__(self, key, value):
        self._parent._load()

        # treat None as delete
        if value is None:
            self._parent._del(self._section, key)
            return

        self._parent._raw.setdefault(self._section, {})
        self._parent._config.setdefault(self._section, {})

        self._parent._raw[self._section][key] = value
        self._parent._config[self._section][key] = value

        self._parent._save()

    def __delitem__(self, key):
        self._parent._del(self._section, key)

    def keys(self):
        return self._parent._config[self._section].keys()

    def items(self):
        return self._parent._config[self._section].items()

    def values(self):
        return self._parent._config[self._section].values()


class _KeyProxy:
    """
    Acts as a live proxy to a key inside a section of a YAMLConfig
    """

    def __init__(self, parent: "YAMLConfig", section: str, key: str):
        self._parent = parent
        self._section = section
        self._key = key

    def __getitem__(self, param):
        self._parent._load()
        return self._parent._config[self._section][self._key].get(param)

    def __setitem__(self, param, value):
        self._parent._load()

        self._parent._raw[self._section].setdefault(self._key, {})
        self._parent._config[self._section].setdefault(self._key, {})

        self._parent._raw[self._section][self._key][param] = value
        self._parent._config[self._section][self._key][param] = value

        self._parent._save()

    def __repr__(self):
        return repr(self._parent._config[self._section][self._key])


def translate_dict(raw: dict, translator: FunctionType) -> dict:
    """
    Applies a translation function to a yaml config (nested dictionary)
    """

    result = {}
    for section, values in raw.items():
        if isinstance(values, dict):
            result[section] = {k: translator(v) for k, v in values.items()}
        else:
            result[section] = translator(values)

    return result


class YAMLConfig:
    """
    Base class for a ConfigParser/dict-style configuration manager
    """

    def __init__(
        self,
        path: Path,
        defaults: dict,
        translator: FunctionType | None = None,
        processor: FunctionType | None = None,
        printer: FunctionType | None = None,
    ):
        self._path = path
        self._defaults = defaults
        self._translator = translator or (lambda x: x)
        self._processor = processor
        self._printer = printer or default_printer

        self._raw = None
        self._config = None
        self._load()

    def _load(self):
        """
        Load yaml and dict representations of config from file if it exists, if not create one with default values. Optionally applies given processing function to dictionary to modify its structure, and can take a printer to modify show() behaviour
        """

        if self._path.exists():
            with open(self._path) as f:
                raw = yaml.safe_load(f) or {}
        else:
            raw = self._defaults
            self._save()

        self._raw = raw
        self._config = translate_dict(raw, self._translator)

        # optional transformation
        if self._processor:
            self._config = self._processor(self._config)

    def _save(self):
        """
        Save yaml representation of config to file
        """

        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "w") as f:
            yaml.dump(self._raw, f, sort_keys=False, indent=YAML_INDENT, Dumper=CustomDumper)

    def reset(self):
        """
        Reset yaml and dict representation of config to their default states, and write this to the config file
        """

        with open(self._path, "w") as f:
            yaml.dump(self._defaults, f, sort_keys=False, indent=YAML_INDENT, Dumper=CustomDumper)
        self._load()

    def _set(self, section: str, key: str, **kwargs):
        """
        Set the value of a given config key in a given section
        """

        self._load()
        if section not in self._raw:
            raise ValueError(f"Section '{section}' not found in config.")

        if section not in self._raw:
            self._raw[section] = {}
        if key not in self._raw[section]:
            self._raw[section][key] = {}

        if section not in self._config:
            self._config[section] = {}
        if key not in self._config[section]:
            self._config[section][key] = {}

        for kwarg, val in kwargs.items():
            if val is not None:
                self._raw[section][key][kwarg] = val
                self._config[section][key][kwarg] = val

        self._save()

    def _del(self, section: str, key: str) -> None:
        """
        Deletes a given key in a given section of the config
        """

        self._load()
        if section not in self._raw:
            raise ValueError(f"Section '{section}' not found in config.")
        if key not in self._raw[section]:
            raise ValueError(f"Key '{key}' not found in section '{section}'.")
        del self._config[section][key]
        del self._raw[section][key]
        self._save()

    def get(self, section: str, key: str, fallback=None):
        """
        Get the value of a key in a given section
        """

        self.load()
        return self._(section, {}).get(key, fallback)

    def get_section(self, section: str):
        """
        Get an entire section by its key in the config
        """

        self._load()
        return self._(section, {})

    def _as_dict(self):
        """
        Get entire config as a 2D dictionary
        """

        self._load()
        return self._config

    def show(self):
        """
        Print config to stdout
        """

        print(self._printer(self._config))

    def open(self) -> None:
        """
        Opens the config file in the default text editor
        """

        open_file(self._path)

    def __getitem__(self, section: str):
        self._load()

        if section not in self._config:
            raise KeyError(f"Section '{section}' not found.")

        return _SectionProxy(self, section)
