import warnings
from enum import Enum, auto

import astropy.units as u
from astropy.units import Quantity

from ..configuration.base_config import BASE_CONFIG


class REQUIRED(Enum):
    """
    REQUIRED.LATER = checks are performed later in query process (e.g. atlas username and password are checked only in atlas lightcurve queries.
    REQUIRED.NOW = needed immediately, i.e. query will always fail without it.
    """

    NOW = auto()
    LATER = auto()


# map of necessary arguments and their default values for each query type. Those with config values do not need to be provided by the user
QUERY_ARGUMENTS = {
    # one of 'survey' and 'catalogue' needed
    "vizier": {"survey": REQUIRED.NOW, "radius": BASE_CONFIG._get("query_settings", "query_radius")},
    # ATLAS requires username and password
    "lightcurve": {
        "survey": REQUIRED.NOW,
        "username": REQUIRED.LATER,
        "password": REQUIRED.LATER,
        "radius": BASE_CONFIG._get("query_settings", "query_radius"),
        "split": False,
        "filter": True,
    },
    "image": {"survey": REQUIRED.NOW, "size": BASE_CONFIG._get("query_settings", "image_size"), "overlays": None, "band": REQUIRED.NOW},
    "spectrum": {"survey": REQUIRED.NOW, "radius": BASE_CONFIG._get("query_settings", "query_radius")},
    # correction needs to be deferred as SED queries use data queries under-the-hood
    "sed": {"radius": BASE_CONFIG._get("query_settings", "query_radius"), "defer_correction": True},
    "hrd": {"survey": "gaia", "colour": "BPmag-RPmag", "mag": "Gmag", "defer_correction": True},
    "datatable": {"columns": REQUIRED.NOW, "radius": BASE_CONFIG._get("query_settings", "query_radius"), "defer_correction": True},
}

UNIVERSAL_ARGUMENTS = {"path": None}

# get default unit scale from config
default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")

REQUIRED_UNITS = {"radius": default_unit, "size": default_unit}


def set_quantity(out_args: dict, parameter: str, unit: Quantity) -> dict:
    value = out_args.get(parameter)
    if value is None:
        return out_args

    if not isinstance(value, Quantity):
        out_args[parameter] = value * unit

    return out_args


def get_query_arguments(kind: str, kwargs: dict) -> dict:
    """
    Checks the arguments for a given query kind, filling in from default (i.e. config) values where possible. Returns the updated dict of kwargs.
    """

    defaults = QUERY_ARGUMENTS.get(kind)
    if defaults is None:
        raise ValueError(f"Invalid query kind '{kind}'. Accepted query kinds: {', '.join(QUERY_ARGUMENTS.keys())}.")
    for key, val in UNIVERSAL_ARGUMENTS.items():
        defaults[key] = val

    missing = [f"'{key}'" for key, val in defaults.items() if val is REQUIRED.NOW and key not in kwargs]
    if missing:
        raise ValueError(f"Missing required arguments(s) for {kind} query: {', '.join(missing)}")

    # output arguments
    out_args = dict(kwargs)

    # any exceptions to the normal rules
    # e.g. ATLAS takes no radius parameter, but setting this to REQUIRED.LATER would require checking that a radius was provided for every other survey
    if kind == "lightcurve":
        # ATLAS doesn't take a radius
        if out_args.get("survey") == "atlas":
            if out_args.get("radius"):
                warnings.warn(
                    "ATLAS light curves are provided as forced photometry at an exact position, and hence setting the radius will have no effect."
                )
            # ATLAS doesn't have object IDs
            if out_args.get("split"):
                warnings.warn(
                    "ATLAS light curves are provided as forced photometry, and hence object IDs to not apply and no splitting will be performed."
                )

            # ATLAS does its own proper motion correction
            out_args["defer_correction"] = True

    # iterate through parameters, setting defaults from the config if not marked as REQUIRED.NOW or REQUIRED.LATER
    for key, val in defaults.items():
        if val not in [REQUIRED.NOW, REQUIRED.LATER]:
            out_args.setdefault(key, val)
        else:
            out_args.setdefault(key, None)

    # set astropy units where required
    for key, val in REQUIRED_UNITS.items():
        out_args = set_quantity(out_args, key, val)

    return out_args
