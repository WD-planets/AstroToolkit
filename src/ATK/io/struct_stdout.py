import inspect
import re

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits.hdu import PrimaryHDU
from astropy.time import Time
from astropy.wcs import WCS

from ..structures.data_containers.definitions import Image

# GLOBALS for tracking state
CURRENT_DEPTH = 0
COL_WIDTHS = {0: 0}
OUTPUT = ""

# OPTIONS
SIG_FIGS = 3  # significant figures of array elements
MAX_DISPLAY = 10  # max entries in an array before truncation occurs
METHODS_TO_IGNORE = ["__eq__", "__init__", "__repr__", "__str__"]

# Headers are printed for containers that need to be expanded
CONTAINER_HEADERS = {pd.DataFrame: lambda x: "<pandas.DataFrame>", dict: lambda x: "<dict>", Image: lambda x: f"<{x.__repr__()}>"}

# ------------------
# SPECIAL FORMATTERS
# ------------------


def format_skycoord(coord: SkyCoord) -> str:
    return f"{round(coord.ra.deg, 3)}°, {round(coord.dec.deg, 3)}°"


def format_dict(dct: dict) -> str:
    """
    Format dict recursively into string representation
    """

    str_rep = ""

    # If dict empty, append label + {} and return
    if not dct:
        return "{}"

    # compute padding for keys inside dict
    key_pad = get_dict_pad(dct, "key: ")

    for key, val in dct.items():
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
        return "[]"

    # recursively format each item
    for item in lst:
        str_rep += format_value(item)

    return str_rep


def format_array(array: np.ndarray | pd.Series) -> str:
    """
    Format array into string representation
    """

    array = np.asarray(array)

    # format array
    formatted = []
    for v in array:
        if v == "...":
            formatted.append("...")
        elif isinstance(v, (float, np.floating)):
            formatted.append(f"{v:.{SIG_FIGS}g}")
        else:
            formatted.append(str(v))

    # truncate array values if length beyond MAX_DISPLAY
    if len(formatted) > MAX_DISPLAY:
        half = MAX_DISPLAY // 2
        formatted = formatted[:half] + ["..."] + formatted[-half:]

    str_rep = "[" + ", ".join(formatted) + "]"

    return str_rep


# ----------------
# SPECIAL TYPE MAP
# ----------------

# Maps to special formatting functions
SPECIAL_FORMATTERS = {
    dict: format_dict,
    pd.DataFrame: lambda df: format_dict(dataframe_to_np_dict(df)),
    np.ndarray: format_array,
    pd.Series: format_array,
    SkyCoord: format_skycoord,
    Time: lambda t: f"{t}",
    PrimaryHDU: lambda h: "<PrimaryHDU>",
    WCS: lambda w: "<WCS>",
}

# -------
# UTILITY
# -------


def is_special_type(value: any) -> bool:
    """
    Check if value has a special formatter function
    """

    if type(value) in SPECIAL_FORMATTERS:
        return True


def apply_special_formatter(value: any) -> str:
    """
    Applies a special formatter function to value
    """

    formatter = SPECIAL_FORMATTERS.get(type(value))
    if formatter:
        return formatter(value)


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
    Replaces placeholder columns in the output str using regex. Placeholders ending in 'H' are placed at half the justification for type hints
    """

    depth = int(match.group(1))
    half = match.group(2) == "H"

    pad = 0
    for key, val in COL_WIDTHS.items():
        if key < depth:
            pad += val

    return " " * (pad // 2 if half else pad)


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


def safe_representation(obj: any) -> str:
    """
    Return a one-line representation of an object from __str__ > __repr__ > __name__
    """

    try:
        return CONTAINER_HEADERS[type(obj)](obj)
    except Exception:
        pass

    for getter in (str, repr):
        try:
            return getter(obj).splitlines()[0]
        except Exception:
            pass

    return f"<{obj.__class__.__name__}>"


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
    Dispatches values to formatters to
    """

    global CURRENT_DEPTH, OUTPUT

    if type(value) in CONTAINER_HEADERS:
        line = f"\n{pad_placeholder(CURRENT_DEPTH + 1, True)}{safe_representation(value)}\n"
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

    methods = [name for name, f in inspect.getmembers(cls, inspect.ismethod) if name not in METHODS_TO_IGNORE and not inspect.isbuiltin(f)]

    return "Available Methods: " + ", ".join(f".{m}()" for m in methods)


def pprint_structure(structure: any, show_all_types: bool) -> None:
    """
    Prints a structure's attributes and methods in a human-readable format. Optionally also prints the types of attributes.
    """

    global CURRENT_DEPTH, OUTPUT

    # get structure attrs
    attrs = structure.__dict__

    # get list of types in attributes
    if show_all_types:
        attrs = add_types_to_keys(attrs)

    # calculate base pad
    pad = get_dict_pad(attrs, ".attr: ")

    inherited_attrs, own_attrs = split_instance_attributes(structure)

    for index, attr_group in enumerate([inherited_attrs, own_attrs]):
        if show_all_types:
            attr_group = add_types_to_keys(attr_group)

        for attr, val in attr_group.items():
            line = f".{attr}: ".ljust(pad)

            update_col_widths(len(line), CURRENT_DEPTH)

            OUTPUT += line
            OUTPUT += format_value(val) + "\n"

        if not index:
            OUTPUT += "\n"

    OUTPUT += print_methods(structure)

    re_exp = re.compile(r"<\|LPAD_DEPTH_(\d+)(H?)\|>")
    formatted = re_exp.sub(replace_pad_placeholder, OUTPUT)

    print(formatted)
