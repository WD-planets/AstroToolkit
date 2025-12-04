import importlib
import inspect
import pkgutil


def build_map(root_module: str, function_name: str, **kwargs):
    """
    Builds a map {str:function} by searching through a given folder for modules that include a given prefix or suffix and contain a known, shared function name
    """

    prefix, suffix = kwargs.get("prefix", None), kwargs.get("suffix", None)
    if prefix is None == suffix is None:
        raise ValueError("build_map() requires exactly one of 'prefix', 'suffix'.")

    map = {}

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
            map[name] = getattr(imported_module, function_name)
        else:
            raise ImportError(f"Module {module.name}.py lacks required function '{function_name}'.")

    return map
