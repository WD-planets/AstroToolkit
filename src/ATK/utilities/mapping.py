import importlib
import inspect
import pkgutil


def build_map(module: str, function_name: str, **kwargs):
    """
    Builds a map {str:function} by searching through a given folder for modules that include a given prefix or suffix and contain a known, shared function name
    """

    prefix, suffix = kwargs.get("prefix", None), kwargs.get("suffix", None)
    if prefix is None == suffix is None:
        raise ValueError("build_map() requires exactly one of 'prefix', 'suffix'.")

    map = {}

    for module_info in pkgutil.iter_modules(module.__path__):
        name = module_info.name

        if prefix and not name.startswith(prefix):
            continue
        if suffix and not name.endswith(suffix):
            continue

        module = importlib.import_module(f"{module.__name__}.{name}")

        if prefix:
            name = name[len(prefix) :]
        elif suffix:
            name = name[: -len(suffix)]

        if hasattr(module, function_name):
            map[name] = getattr(module, function_name)
        else:
            raise ImportError(f"Module {name} lacks required function '{name}'.")

    return map


def build_structure_map():
    module = importlib.import_module("ATK.structures.query_results.definitions")
    suffix = "Struct"

    map = {}
    for name, obj in inspect.getmembers(module, inspect.isclass):
        if obj.__module__ == module.__name__ and name.endswith(suffix):
            key = name[: -len(suffix)]
            map[key.lower()] = obj

    return map
