import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.units import Quantity
from scipy import stats

np.seterr(divide="ignore")
warnings.simplefilter("ignore", category=RuntimeWarning)


def bin_by_size(x: np.ndarray, ys: list[np.ndarray], y_errs: list[np.ndarray], size: float):
    t_min = np.nanmin(x)
    t_max = np.nanmax(x)
    edges = np.arange(t_min, t_max + size, size)

    return do_binning(x, ys, y_errs, edges)


def do_binning(x: np.ndarray, ys: list[np.ndarray], y_errs: list[np.ndarray], bins: int | np.ndarray):
    x_bins, out_ys, out_y_errs = None, [], []

    for y, err in zip(ys, y_errs):
        save_errors = True
        if err is None:
            err = np.ones(x.shape)
            save_errors = False

        weights = 1 / np.power(err, 2)
        # numerator
        numerator, bin_edges, _ = stats.binned_statistic(x, y * weights, statistic="sum", bins=bins)
        # denominator
        denominator, _, _ = stats.binned_statistic(x, weights, statistic="sum", bins=bin_edges)

        w_mean = numerator / denominator
        w_mean_err = np.sqrt(1.0 / denominator)
        mid_bins = (bin_edges[1:] + bin_edges[:-1]) / 2
        df = pd.DataFrame({"w_mean": w_mean, "w_mean_err": w_mean_err, "mid_bins": mid_bins})
        df = df.dropna()

        w_mean, w_mean_err, x_bins = df.T.to_numpy()

        out_ys.append(w_mean)
        if save_errors:
            out_y_errs.append(w_mean_err)

    return x_bins, out_ys, out_y_errs


def bin_nd(
    x: np.ndarray,
    ys: list[np.ndarray],
    errs: list[np.ndarray] | None = None,
    bins: int | None = None,
    size: Quantity | float | None = None,
):
    if bins is None == size is None:
        raise ValueError("Exactly one of 'bins', 'size' must be provided.")

    while len(errs) < len(ys):
        errs.append(None)

    if size:
        if isinstance(size, Quantity):
            size = size.value
        out_x, out_ys, out_y_errs = bin_by_size(x, ys, errs, size)
    else:
        out_x, out_ys, out_y_errs = do_binning(x, ys, errs, bins)

    return out_x, out_ys, out_y_errs
