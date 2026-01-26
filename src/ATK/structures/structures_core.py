import copy
from dataclasses import dataclass

import numpy
import pandas
from astropy.io.fits import BinTableHDU
from astropy.units import Quantity, Unit

from .methods.lightcurve.phasefold import fold_lc
from .methods.lightcurve.powspec import gen_powspec

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
class BaseContainer:
    _target_key: str | None = None

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

    def __repr__(self):
        return f"<{self.survey} {type(self).__name__}>"

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
        from ..io.structure_io import struct_to_dataframe

        return struct_to_dataframe(self)

    @classmethod
    def from_dataframe(cls, data: pandas.DataFrame, **kwargs: any):
        from ..io.structure_io import struct_from_dataframe

        return struct_from_dataframe(cls, data, **kwargs)

    def to_hdu(self) -> BinTableHDU:
        from ..io.structure_io import struct_to_hdu

        return struct_to_hdu(self)
