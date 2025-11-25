import importlib

from astropy.coordinates import SkyCoord

from .configuration.base_config import BASE_CONFIG
from .queries.arguments import get_query_arguments
from .utilities.defaults import RETURNS
from .utilities.mapping import build_map
from .utilities.targeting import prepare_search


def query(kind: str, target: int | SkyCoord, **kwargs):
    module = importlib.import_module(f"ATK.queries.{kind}")
    query_map = build_map(module, "query", suffix="_query")

    # get necessary parameters from config if not given
    additional_arguments = get_query_arguments(kind, kwargs)

    # get specific query function (for given survey if multiple are available)
    query_function = query_map[kwargs.get("survey")] if len(query_map) > 1 else list(query_map.values())[0]

    # get search position and structure
    search_pos, structure = prepare_search(target=target, query_kind=kind, **additional_arguments)

    # add source to kwargs (needed e.g. in Vizier queries to Gaia by source)
    if isinstance(target, int):
        additional_arguments["source"] = target

    # perform query
    query_result = query_function(search_pos, **additional_arguments)

    if query_result is RETURNS.EXCEPTION:
        structure.data = None
        structure.exception = True
    elif query_result is RETURNS.NULL:
        structure.data = None
    else:
        structure.data = query_result

    return structure
