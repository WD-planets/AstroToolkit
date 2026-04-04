import numpy as np
from astropy.stats import sigma_clip


def do_sigma_clipping(y: np.ndarray, arrs: list[np.ndarray], sigma: float, sigma_lower: float = None, sigma_upper: float = None):
    mask = sigma_clip(data=y, sigma=sigma, sigma_lower=sigma_lower, sigma_upper=sigma_upper, masked=True).mask
    keep = ~mask

    y = y[keep]
    for i, arr in enumerate(arrs):
        arrs[i] = arr[keep]

    # Need to safeguard against no data remaining?

    return y, arrs
