import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.units import Quantity
from scipy import stats


def _strip_unit(arr):
    if isinstance(arr, Quantity):
        return arr.value, arr.unit
    return arr, None


def _align_to_unit(input, target_unit, input_name: str, ref_name: str):
    _, in_unit = _strip_unit(input)

    if target_unit is None and in_unit is None:
        return input
    if target_unit is not None and in_unit is None:
        return input
    if target_unit is None and in_unit is not None:
        raise u.UnitConversionError(f"'{input_name}' is a Quantity (unit: '{in_unit}') but '{ref_name}' has no unit.")

    try:
        return input.to(target_unit).value
    except u.UnitConversionError:
        raise u.UnitConversionError(f"Cannot convert '{input_name}' (unit: '{in_unit}') to '{ref_name}' unit '{target_unit}'.")
