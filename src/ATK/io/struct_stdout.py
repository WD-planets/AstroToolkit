import inspect
import re
from enum import Enum

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import BinTableHDU, ImageHDU, PrimaryHDU
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity, UnrecognizedUnit
from astropy.wcs import WCS
from bokeh.layouts import GridBox
from bokeh.models import Column, Row
from bokeh.plotting import figure

from ..configuration.base_config import BASE_CONFIG
from ..structures.Target import Target
from ..utilities.mapping import build_structure_map

# this should be left to False, kwarg 'debug' can be used to set it locally
DEBUG = False
SHOW_ALL = False
SHOW_TYPES = False

STRUCTURE_MAP = build_structure_map()
NO_COMMA_SEP = list(STRUCTURE_MAP.values()) + [pd.DataFrame]

# GLOBALS for tracking state
CURRENT_DEPTH = 0
COL_WIDTHS = {0: 0}
OUTPUT = ""
SEEN_IDS = set()

# OPTIONS
ROUND = 3  # decimal places of float elements
MAX_DISPLAY = 4  # max entries in an array before truncation occurs

# Headers are printed for containers that need to be expanded (ATK containers are handled separately using MAP
CONTAINER_HEADERS = {pd.DataFrame: lambda x: "<pandas.DataFrame>", dict: lambda x: "<dict>", Table: "<astropy.Table>"}

# ------------------
# SPECIAL FORMATTERS
# ------------------

UNITS = {u.arcsec: "″", u.arcmin: "′", u.deg: ["°", " deg"]}
unit_format = BASE_CONFIG._get("global_settings", "unit_format")
if unit_format == "symbol":
    UNITS = {key: val if not isinstance(val, list) else val[0] for key, val in UNITS.items()}
elif unit_format == "text":
    UNITS = {key: val[1] for key, val in UNITS.items() if isinstance(val, list)}
else:
    raise ValueError(f"Unexpected config entry for unit_format '{unit_format}'.")


def format_target(target: Target) -> str:
    """
    Format ATK Target into string representation
    """

    if target.identifier:
        str_rep = f"{target.identifier} | "
    else:
        str_rep = ""

    str_rep += format_skycoord(target.initial_coords)

    if target.radius is not None:
        str_rep = f"{str_rep[:-1]}, {target.radius.value}{UNITS[target.radius.unit]})"

    return str_rep


def format_skycoord(coord: SkyCoord) -> str:
    """
    Format SkyCoord into string representation
    """

    str_rep = f"{round(coord.ra.deg, ROUND)}{UNITS[u.deg]} {round(coord.dec.deg, ROUND)}{UNITS[u.deg]}"

    if hasattr(coord, "frame") and hasattr(coord, "obstime"):
        str_rep += f" ({coord.frame.name}, {coord.obstime.fits})"
    elif hasattr(coord, "frame"):
        str_rep += f" ({coord.frame.name})"
    elif hasattr(coord, "obstime"):
        str_rep += f" ({coord.obstime.fits})"

    return str_rep


def format_dict(dct: dict, show_types_override: bool | None = None, origin_type=dict) -> str:
    """
    Format dict recursively into string representation
    """

    str_rep = ""

    # If dict empty, append label + {} and return
    if not dct:
        return "{}"

    if (SHOW_TYPES or show_types_override is True) and (show_types_override is None or show_types_override is True):
        dct = add_types_to_keys(dct)

    # compute padding for keys inside dict
    key_pad = get_dict_pad(dct, "key: ")

    for key, val in dct.items():
        if key.startswith("_") and not DEBUG and origin_type is not pd.DataFrame:
            continue

        if val is None:
            continue

        str_rep += pad_placeholder(CURRENT_DEPTH)

        line = f"{key}: ".ljust(key_pad)
        update_col_widths(len(line), CURRENT_DEPTH)
        str_rep += line

        # Recurse for each value
        str_rep += format_value(val) + "\n"

    return str_rep


def format_list(lst: list) -> str:
    """
    Formats lists recursively into string representation. Arrays should be used for simple objects.
    """

    str_rep = ""

    if not lst:
        return "<empty list>\n"

    # recursively format each item
    for i, item in enumerate(lst):
        if i == MAX_DISPLAY and not SHOW_ALL:
            str_rep += format_value(f"+{len(lst) - MAX_DISPLAY} more ...")
            return str_rep

        str_rep += format_value(item)
        if type(item) not in NO_COMMA_SEP and i != len(lst) - 1:
            str_rep += ", "

    return str_rep


def format_array(array: np.ndarray | pd.Series) -> str:
    """
    Format array into string representation
    """

    array = np.asarray(array)
    str_rep = "[]"

    # format array
    formatted = []
    for v in array:
        if isinstance(v, np.ndarray):
            formatted.append(format_value(v))
        elif v == "...":
            formatted.append("...")
        elif isinstance(v, (float, np.floating)):
            formatted.append(f"{round(v, ROUND)}")
        else:
            formatted.append(str(v))

    # truncate array values if length beyond MAX_DISPLAY
    if len(formatted) > MAX_DISPLAY:
        half = MAX_DISPLAY // 2
        formatted = formatted[:half] + ["..."] + formatted[-half:]

    str_rep = str_rep[0] + ", ".join(formatted) + str_rep[-1]

    return str_rep


def format_quantity(val: Quantity) -> str:
    """
    Formats astropy Quantities or array-like Quantities into string representation
    """

    if val.shape:
        val_arr = format_value(val.value)
    else:
        val_arr = f"{round(val.value, ROUND)}"

    if isinstance(val.unit, UnrecognizedUnit):
        return f"{val_arr} {val.unit.to_string('unicode')}"

    if val.unit in UNITS:
        str_rep = f"{val_arr}{' ' if val.shape else ''}{UNITS[val.unit]}"
    elif unit_format == "text":
        str_rep = f"{val_arr} {val.unit.to_string()}"
    elif unit_format == "symbol":
        str_rep = f"{val_arr} {val.unit.to_string('unicode')}"

    return str_rep


def format_time(time: Time) -> str:
    return f"{time.fits}"


# ----------------
# SPECIAL TYPE MAP
# ----------------

# Maps to special formatting functions
SPECIAL_FORMATTERS = {
    dict: format_dict,
    pd.DataFrame: lambda df: format_dict(dataframe_to_np_dict(df), False, pd.DataFrame),
    Table: lambda tbl: format_dict(table_to_dict(tbl), False, Table),
    np.ndarray: format_array,
    pd.Series: format_array,
    Target: format_target,
    SkyCoord: format_skycoord,
    Time: format_time,
    PrimaryHDU: lambda x: "<PrimaryHDU>",
    ImageHDU: lambda x: "<ImageHDU>",
    BinTableHDU: lambda x: "<BinTableHDU>",
    WCS: lambda x: "<WCS>",
    figure: lambda x: "<Bokeh Figure>",
    Row: lambda x: "<Bokeh Figure>",
    Column: lambda x: "<Bokeh Figure>",
    GridBox: lambda x: "<Bokeh Figure>",
    Quantity: format_quantity,
}

# -------
# UTILITY
# -------


def is_special_type(value: any) -> bool:
    """
    Check if value has a special formatter function, or if it is an enum
    """

    if type(value) in SPECIAL_FORMATTERS:
        return True

    if isinstance(value, Enum):
        return True


def apply_special_formatter(value: any) -> str:
    """
    Applies a special-case formatter function to value
    """

    formatter = SPECIAL_FORMATTERS.get(type(value))
    if formatter:
        return formatter(value)

    # Enum fallback to avoid recursion error + print name of enum value
    if isinstance(value, Enum):
        return value.name


def get_dict_pad(dictionary: dict, key_example: str) -> int:
    """
    Calculated the padding needed for a dict object, using a given representation to add additional padding for any non-key characters
    """

    return max(len(attr) for attr in dictionary) + len(key_example) - sum(map(str.isalpha, key_example))


def update_col_widths(pad: int, depth: int) -> None:
    """
    Updates the width of the current column to the maximum size of a key in that column (i.e. at the current recursion depth)
    """

    global COL_WIDTHS

    if pad > COL_WIDTHS.get(depth, -1):
        COL_WIDTHS[depth] = pad


def pad_placeholder(depth: int, half: bool = False) -> str:
    """
    Generates a placeholder for padding at the current recursion level, this is later replaced by white space to left-justify keys
    """

    return f"<|LPAD_DEPTH_{depth}{'H' if half else ''}|>"


def replace_pad_placeholder(match: re.Match) -> str:
    """
    Replaces placeholder columns in the output str using regex. Placeholders ending in 'H' are placed at half the justification (for container types)
    """

    depth = int(match.group(1))
    half = match.group(2) == "H"

    pad = 0
    last_pad = 0
    for key, val in COL_WIDTHS.items():
        if key < depth:
            pad += val
            last_pad = val
    if half:
        final_pad = pad - last_pad // 2
    else:
        final_pad = pad

    return " " * final_pad


def is_expandable(val: any) -> bool:
    """
    Check if an object's attributes should be expanded, special types excluded
    """

    return isinstance(val, dict) or hasattr(val, "__dict__")


def dataframe_to_np_dict(df: pd.DataFrame) -> dict:
    """
    Converts a dataframe to a dict with all columns converted to numpy arrays
    """

    dct = df.to_dict(orient="list")
    for key, val in dct.items():
        if isinstance(val, list):
            dct[key] = np.asarray(val)

    return dct


def table_to_dict(tbl: Table) -> dict:
    """
    Converts an astropy Table into a dict with all columns converted to numpy arrays
    """

    out = {}

    out = {}

    for name in tbl.colnames:
        col = tbl[name]

        # Get plain array (handles MaskedColumn safely)
        arr = np.array(col)

        # If masked, replace mask in a dtype-safe way
        if hasattr(col, "mask") and np.any(col.mask):
            mask = col.mask

            if np.issubdtype(col.dtype, np.number):
                # numeric → safe float representation
                arr = arr.astype(float)
                arr[mask] = np.nan
            else:
                # strings / objects → keep as object, fill None
                arr = arr.astype(object)
                arr[mask] = None

        # Attach units if present
        unit = getattr(col, "unit", None)
        if unit is not None:
            out[name] = Quantity(arr, unit=unit)
        else:
            out[name] = arr

    return out


def safe_representation(obj: any) -> str:
    """
    Return a one-line representation of an object. Unless overriden, this is (type) for builtins and (package.type) for others
    """

    try:
        return CONTAINER_HEADERS[type(obj)](obj)
    except Exception:
        pass

    module_base = type(obj).__module__.split(".")[0]
    typ_name = type(obj).__name__

    if module_base == "builtins":
        return f"({typ_name})"

    return f"({module_base}.{typ_name})"


def add_types_to_keys(dct: dict):
    """
    Adds the type of each value in a dict to its corresponding key (example: 1 -> example (int): 1)
    """

    typed_dict = {}
    for key, val in dct.items():
        rep = safe_representation(val)
        if type(val) not in CONTAINER_HEADERS:
            typed_dict[f"{key} {rep}"] = val
        else:
            typed_dict[key] = val

    return typed_dict


# ----------
# DISPATCHER
# ----------


def format_value(value: any) -> str:
    """
    Dispatches values to formatters
    """

    global CURRENT_DEPTH, OUTPUT, SEEN_IDS

    # recursion guard
    obj_id = id(value)
    if obj_id in SEEN_IDS:
        return f"<recursion:{type(value).__name__}>"
    SEEN_IDS.add(obj_id)

    # 3rd party expandable objects
    if type(value) in CONTAINER_HEADERS:
        line = f"\n{pad_placeholder(CURRENT_DEPTH + 1, True)}{safe_representation(value)}\n"
    # ATK objects
    elif type(value) in STRUCTURE_MAP.values():
        line = f"\n{pad_placeholder(CURRENT_DEPTH + 1, True)}{value.__repr__()}\n"
    # everything else
    else:
        line = ""

    try:
        CURRENT_DEPTH += 1
        if is_special_type(value):
            line += f"{apply_special_formatter(value)}"
        elif isinstance(value, list):
            line += f"{format_list(value)}"
        elif is_expandable(value):
            line += format_dict(value.__dict__)
        else:
            line += f"{value}"

        return line

    finally:
        SEEN_IDS.discard(obj_id)
        CURRENT_DEPTH -= 1


# ----------------
# CLASS INSPECTION
# ----------------


def split_instance_attributes(obj: any) -> tuple[dict]:
    """
    Splits the attributes of a class into its own and those of its parent
    """

    cls = type(obj)
    inherited, own = {}, {}

    base_attrs = set()
    for base in cls.__mro__[1:]:
        base_attrs |= set(getattr(base, "__annotations__", {}).keys())

    for key, val in vars(obj).items():
        if key in base_attrs:
            inherited[key] = val
        else:
            own[key] = val

    return inherited, own


# ----
# MAIN
# ----


def print_methods(cls: any) -> str:
    """
    Prints available methods of an object, excluding
    """

    methods = [name for name, f in inspect.getmembers(cls, inspect.ismethod) if not inspect.isbuiltin(f) and not name.startswith("_")]

    return "\nAvailable Methods: " + ", ".join(f".{m}()" for m in methods)


def pprint_structure(structure: any, show_types: bool, **kwargs) -> None:
    """
    Prints a structure's attributes and methods in a human-readable format. Optionally also prints the types of attributes.
    """

    global CURRENT_DEPTH, OUTPUT, COL_WIDTHS, DEBUG, SHOW_ALL, SHOW_TYPES

    DEBUG = kwargs.get("debug", False)
    SHOW_ALL = kwargs.get("show_all", False)
    SHOW_TYPES = show_types

    # get structure attrs
    attrs = structure.__dict__

    # get list of types in attributes
    if show_types:
        attrs = add_types_to_keys(attrs)

    # calculate base pad
    pad = get_dict_pad(attrs, ".attr: ")

    print(structure.__repr__())

    # split attributes into inherited + uninherited
    inherited_attrs, own_attrs = split_instance_attributes(structure)

    for index, attr_group in enumerate([inherited_attrs, own_attrs]):
        # add type strings
        if show_types or DEBUG:
            attr_group = add_types_to_keys(attr_group)

        # iterate through attributes in group if requested
        for attr, val in attr_group.items():
            if not DEBUG and (val is None or attr.startswith("_")):
                continue

            line = f".{attr}: ".ljust(pad)

            update_col_widths(len(line), CURRENT_DEPTH)

            # useful for debugging, e.g. showing where a recursion depth error was encountered
            if DEBUG:
                print(attr, val)

            OUTPUT += line
            OUTPUT += format_value(val) + "\n"

        if not index:
            OUTPUT += "\n"

    OUTPUT += print_methods(structure)

    # replace placeholder strings with whitespace
    re_exp = re.compile(r"<\|LPAD_DEPTH_(\d+)(H?)\|>")
    formatted = re_exp.sub(replace_pad_placeholder, OUTPUT)

    # reset global variables
    OUTPUT = ""
    CURRENT_DEPTH = 0
    COL_WIDTHS = {0: 0}

    formatted += "\n"

    print(formatted)
