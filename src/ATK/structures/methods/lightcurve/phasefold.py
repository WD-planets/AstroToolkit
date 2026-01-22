from __future__ import annotations

from typing import TYPE_CHECKING

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.definitions import Lightcurve


def fold_lc(struct: object, lcs: list[Lightcurve], min: float, max: float, samples: int):
    freq = do_ls(lcs, min, max, samples)
