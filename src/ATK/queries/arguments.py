from enum import Enum, auto

from ..configuration.base_config import BASE_CONFIG


class REQUIRED(Enum):
    """
    REQUIRED.LATER = checks are performed later in query process (e.g. atlas username and password are checked only in atlas lightcurve queries.
    REQUIRED.NOW = needed immediately, i.e. query will always fail without it.
    """

    NOW = auto()
    LATER = auto()


# map of necessary arguments and their default values for each query type
QUERY_ARGUMENTS = {
    "vizier": {"alias": REQUIRED.LATER, "catalogue": REQUIRED.LATER, "radius": BASE_CONFIG.get("query_settings", "query_radius")},
    "lightcurve": {
        "survey": REQUIRED.LATER,
        "username": REQUIRED.LATER,
        "password": REQUIRED.LATER,
        "radius": BASE_CONFIG.get("query_settings", "query_radius"),
    },
    "image": {
        "survey": REQUIRED.NOW,
        "size": BASE_CONFIG.get("query_settings", "image_size"),
        "overlays": REQUIRED.LATER,
        "band": REQUIRED.NOW,
    },
    "spectrum": {"survey": REQUIRED.NOW, "radius": BASE_CONFIG.get("query_settings", "query_radius")},
    "sed": {"radius": BASE_CONFIG.get("query_settings", "query_radius")},
}


def get_query_arguments(kind: str, kwargs: dict) -> dict:
    """
    Checks the arguments for a given query kind, filling in from default (i.e. config) values where possible. Returns the updated dict of kwargs.
    """

    defaults = QUERY_ARGUMENTS[kind]

    missing = [f"'{key}'" for key, val in defaults.items() if val is REQUIRED.NOW and key not in kwargs]
    if missing:
        raise ValueError(f"Missing required arguments(s) for {kind} query: {', '.join(missing)}")

    # output arguments
    out_args = dict(kwargs)

    # iterate through defaults, setting defaults from the config if not marked as REQUIRED.NOW or REQUIRED.LATER
    for key, val in defaults.items():
        if val not in [REQUIRED.NOW, REQUIRED.LATER]:
            out_args.setdefault(key, val)
        else:
            out_args.setdefault(key, None)

    return out_args
