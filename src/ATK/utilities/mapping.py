import importlib
import inspect
import pkgutil
from enum import EnumType
from types import ModuleType

from .defaults import QUERY_KINDS


def get_query_result_map():
    """
    Responds a map {str:class} from query kinds -> query result objects (QueryResult or PlottableQueryResult)
    """

    from ..structures.definitions import PlottableQueryResult, QueryResult

    NON_PLOTTABLE = ["vizier"]

    return {kind: QueryResult if kind in NON_PLOTTABLE else PlottableQueryResult for kind in QUERY_KINDS}


def build_map(root_module: ModuleType, function_name: str, **kwargs):
    """
    Builds a map {str:function} by searching through a given folder for modules that include a given prefix or suffix and contain a known, shared function name
    """

    prefix, suffix = kwargs.get("prefix", None), kwargs.get("suffix", None)
    if prefix is None == suffix is None:
        raise ValueError("build_map() requires exactly one of 'prefix', 'suffix'.")

    func_map = {}

    modules = pkgutil.iter_modules(root_module.__path__)

    for module in modules:
        name = module.name

        if prefix and not name.startswith(prefix):
            continue
        if suffix and not name.endswith(suffix):
            continue

        imported_module = importlib.import_module(f"{root_module.__name__}.{name}")

        if prefix:
            name = name[len(prefix) :]
        elif suffix:
            name = name[: -len(suffix)]

        if hasattr(imported_module, function_name):
            func_map[name] = getattr(imported_module, function_name)
        else:
            raise ImportError(f"Module {module.name}.py lacks required function '{function_name}'.")

    return func_map


def build_structure_map():
    """
    Creates a map of ATK structure definitions, exluding enums
    """

    module = importlib.import_module("ATK.structures.definitions")

    struct_map = {}

    for name, obj in inspect.getmembers(module, inspect.isclass):
        if obj.__module__ == module.__name__ and not isinstance(obj, EnumType):
            struct_map[name] = obj

    return struct_map
