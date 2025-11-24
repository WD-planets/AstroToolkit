import importlib
import inspect
import pkgutil


def build_map(module: str, function_name: str, **kwargs):
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

        if hasattr(module, function_name):
            map[name] = getattr(module, function_name)
        else:
            raise ImportError(f"Module {name} lacks required function '{name}'.")
