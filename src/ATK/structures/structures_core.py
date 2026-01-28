import copy
from dataclasses import dataclass

import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU
from astropy.units import Quantity, Unit

from ..io.structure_io import COLUMN_TYPES, get_cols
from .methods.lightcurve.phasefold import fold_lc
from .methods.lightcurve.powspec import gen_powspec
from .Target import Target

# data methods that required a collection of containers to be processed
GROUP_METHODS = {"fold": fold_lc, "pspec": gen_powspec}

# whether to combine data structures into combined plots
PLOT_METHODS = {
    "image": "individual",
    "lightcurve": "combined",
    "spectrum": "individual",
    "sed": "individual",
    "hrd": "combined",
    "powspec": "individual",
    "phasefold": "combined",
}

# whether to split plots by target
SPLIT_BY_TARGET = {"image": True, "lightcurve": True, "spectrum": True, "sed": True, "hrd": False, "powspec": True, "phasefold": True}

# type hint for arrays of astropy Quantities
QuantityArray = numpy.ndarray[Quantity]


def manage_inplace(structure: any, inplace: bool):
    if inplace:
        return structure
    else:
        if hasattr(structure, "figure") and structure.figure:
            figure = structure.figure
            structure.figure = None
            struct_copy = copy.deepcopy(structure)
            structure.figure = figure

            return struct_copy
        return copy.deepcopy(structure)


@dataclass
class Container:
    _target_key: str | None = None

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

    def __repr__(self):
        survey_str = f"{self.survey} " if self.survey else ""

        return f"<{survey_str}{type(self).__name__}>"

    def __str__(self):
        return self.__repr__()

    def _get_attr_value(self, attr: str) -> numpy.ndarray:
        val = getattr(self, attr)
        if val is None:
            raise ValueError(f"Empty or invalid {type(self).__name__} attribute '{attr}'.")

        if isinstance(val, Quantity):
            return val.value
        return val

    def _get_attr_unit(self, attr: str) -> Unit:
        val = getattr(self, attr)
        if val is None:
            raise ValueError(f"Empty or invalid {type(self).__name__} attribute '{attr}'.")

        if isinstance(val, Quantity):
            return val.unit
        return None

    def _get_cols(self):
        from ..io.structure_io import get_cols

        return get_cols(self)

    def to_dataframe(self) -> pandas.DataFrame:
        return struct_to_dataframe(self)

    @classmethod
    def from_dataframe(cls, target: Target | int | SkyCoord, data: pandas.DataFrame, **kwargs: any):
        return struct_from_dataframe(cls, target, data, **kwargs)

    def to_hdu(self) -> BinTableHDU:
        from ..io.structure_io import struct_to_hdu

        return struct_to_hdu(self)


def struct_from_dataframe(ctnr: any, target: Target | int | SkyCoord, data: pandas.DataFrame, **kwargs: dict) -> any:
    from ..queries.query_core import setup_targeting

    ctnr_cols = get_cols(ctnr)

    relevant_data = {}
    for col in data.columns.values.tolist():
        if hasattr(ctnr, col) and col in ctnr_cols:
            relevant_data[col] = data[col].to_numpy()

    for arg, val in kwargs.items():
        if hasattr(ctnr, arg) and arg not in ctnr_cols:
            relevant_data[arg] = val

    ctnr = ctnr(**relevant_data)

    targets = setup_targeting(target)
    if len(targets) > 1:
        raise ValueError("Only one target may be provided per data container.")

    ctnr._target_key = targets[0]._key

    return ctnr


def struct_to_dataframe(structure: any) -> pandas.DataFrame:
    """
    Combines the array-like attributes of a data structure into a single pandas DataFrame
    """

    cols = get_cols(structure)

    data = {}
    for col in cols:
        val = getattr(structure, col)
        if val is None:
            continue

        if not isinstance(val, COLUMN_TYPES):
            val = [val]
        data[col] = val

    return pandas.DataFrame.from_dict(data)
