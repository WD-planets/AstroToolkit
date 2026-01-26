import importlib
import inspect
import pkgutil
from enum import EnumType
from types import ModuleType

from ..structures.DataSet import DataSet
from ..structures.structures_core import BaseContainer


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

    struct_map = {}

    package = importlib.import_module("ATK.structures")

    for _, module_name, _ in pkgutil.iter_modules(package.__path__):
        module = importlib.import_module(f"{package.__name__}.{module_name}")

        for name, obj in inspect.getmembers(module, inspect.isclass):
            # Ensure class is defined in *this* module
            if obj.__module__ != module.__name__:
                continue

            if isinstance(obj, EnumType):
                continue

            if not issubclass(obj, (BaseContainer, DataSet)):
                continue

            struct_map[name] = obj

    return struct_map
