import importlib

from astropy.coordinates import SkyCoord

from .configuration.base_config import BASE_CONFIG
from .utilities.mapping import build_map


def query(kind: str, target: int | SkyCoord, **kwargs):
    module = importlib.import_module(f"ATK.queries.{kind}")
    query_map = build_map(module, "query", suffix="_query")
    query_function = query_map[kind]

    data = query_function(target)
