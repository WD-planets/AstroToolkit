import numpy as np
from astropy.units import Quantity

from ...utilities.units import _align_to_unit, _strip_unit


def crop_nd(x: np.ndarray, ys: list[np.ndarray], lower_lim: float | None, upper_lim: float | None):
    if lower_lim is None and upper_lim is None:
        raise ValueError("Atleast one of 'lower_lim', 'upper_lim' must be provided.")

    mask = np.ones(x.shape, dtype=bool)

    x, x_unit = _strip_unit(x)
    lower_lim = _align_to_unit(lower_lim, x_unit, "lower_lim", "x_arr")
    upper_lim = _align_to_unit(upper_lim, x_unit, "upper_lim", "x_arr")

    y_units = []
    for i, y in enumerate(ys):
        if isinstance(y, Quantity):
            y_units.append(y.unit)
            ys[i] = y.value
        else:
            y_units.append(None)

    if lower_lim is not None:
        mask &= x >= lower_lim
    if upper_lim is not None:
        mask &= x <= upper_lim

    out_x = x[mask]
    out_ys = [a[mask] for a in ys]

    if x_unit:
        out_x = out_x * x_unit

    for index, (y_unit, out_y) in enumerate(zip(y_units, out_ys)):
        if y_unit is not None:
            out_ys[index] = out_y * y_unit

    return out_x, out_ys
