import numpy as np


def crop_nd(x: np.ndarray, ys: list[np.ndarray], lower_lim: float | None, upper_lim: float | None):
    if not lower_lim and not upper_lim:
        raise ValueError("Atleast one of 'lower_lim', 'upper_lim' must be provided.")

    mask = np.ones_like(x, dtype=bool)

    if lower_lim is not None:
        mask &= x >= lower_lim
    if upper_lim is not None:
        mask &= x <= upper_lim

    out_x = x[mask]
    out_ys = [a[mask] for a in ys]

    return out_x, out_ys
