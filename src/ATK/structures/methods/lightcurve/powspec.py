from __future__ import annotations

from typing import TYPE_CHECKING

from astropy import units as u

from .timeseries_core import do_ls

if TYPE_CHECKING:
    from ....structures.definitions import Lightcurve, Powspec


def gen_powspec(struct: object, lcs: list[Lightcurve], min: float, max: float, samples: int) -> Powspec:
    from ....structures.definitions import Powspec

    struct.kind = "powspec"

    freq, power, fopt = do_ls(lcs, min, max, samples)
    period = (1 / fopt).to(u.day)
    bands = list(set([lc.band for lc in lcs]))

    pspec = Powspec(survey=struct.survey, band=", ".join(bands), frequency=freq, power=power, fopt=fopt, popt=period, _target_key=lcs[0]._target_key)

    return pspec
