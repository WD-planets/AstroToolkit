from enum import Enum, auto

from ..configuration.base_config import BASE_CONFIG


class REQUIRED(Enum):
    NOW = auto()
    LATER = auto()


# map of necessary arguments and their default values for each query type
QUERY_ARGUMENTS = {
    "vizier": {
        "survey": REQUIRED.LATER,
        "catalogue": REQUIRED.LATER,
        "radius": BASE_CONFIG.get("query_settings", "query_radius"),
    },
    "lightcurve": {
        "survey": REQUIRED.LATER,
        "username": REQUIRED.LATER,
        "password": REQUIRED.LATER,
        "radius": BASE_CONFIG.get("query_settings", "query_radius"),
    },
    "image": {
        "survey": REQUIRED.NOW,
        "size": BASE_CONFIG.get("query_settings", "image_size"),
        "overlays": BASE_CONFIG.get("query_settings", "image_overlays"),
        "band": BASE_CONFIG.get("query_settings", "image_band"),
    },
    "spectrum": {"survey": REQUIRED.NOW, "radius": BASE_CONFIG.get("query_settings", "query_radius")},
    "sed": {"radius": BASE_CONFIG.get("query_settings", "query_radius")},
}


def get_query_arguments(kind: str, kwargs):
    defaults = QUERY_ARGUMENTS[kind]

    missing = [key for key, val in defaults.items() if val is REQUIRED.NOW and key not in kwargs]
    if missing:
        raise ValueError(f"Missing required arguments(s) for '{kind}' query: {'. '.join(missing)}")

    out_args = dict(kwargs)
    for key, val in defaults.items():
        if val not in [REQUIRED.NOW, REQUIRED.LATER]:
            out_args.setdefault(key, val)
        else:
            out_args.setdefault(key, None)

    return out_args
