from pathlib import Path
from types import FunctionType

import yaml


def translate_dict(raw: dict, translator: FunctionType) -> dict:
    result = {}
    for section, values in raw.items():
        if isinstance(values, dict):
            result[section] = {k: translator(v) for k, v in values.items()}
        else:
            result[section] = translator(values)

    return result


class YAMLConfig:
    def __init__(self, path: Path, defaults: dict, translator: FunctionType | None = None, processor: FunctionType | None = None):
        self._path = path
        self._defaults = defaults
        self._translator = translator or (lambda x: x)
        self._processor = processor

        self._raw = None
        self._config = None
        self.load()

    def load(self):
        if self._path.exists():
            with open(self._path) as f:
                raw = yaml.safe_load(f) or {}
        else:
            raw = self._defaults

        self._raw = raw
        self._config = translate_dict(raw, self._translator)

        # optional transformation
        if self._processor:
            self._config = self._processor(self._config)

    def save(self):
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "w") as f:
            yaml.dump(self._raw, f, sort_keys=False)

    def reset(self):
        with open(self._path, "w") as f:
            yaml.dump(self._defaults, f, sort_keys=False)
        self.load()

    def set(self, section: str, key: str, value):
        if section not in self._raw:
            raise ValueError(f"Unexpected section '{section}'.")
        self._raw[section][key] = value
        self._config[section][key] = value

    def get(self, section: str, key: str, fallback=None):
        self.load()
        return self._config.get(section, {}).get(key, fallback)

    def get_section(self, section: str):
        return self._config.get(section, {})

    def as_dict(self):
        return self._config

    def show(self):
        for section, values in self._config.items():
            print(f"[{section}]")
            if isinstance(values, dict):
                for key, value in values.items():
                    print(f"{key} = {value}")
            else:
                print(values)
            print()
