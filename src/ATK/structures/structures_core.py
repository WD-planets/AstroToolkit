import copy
from dataclasses import dataclass, field
from typing import Self

import numpy
from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU
from astropy.table import Table
from astropy.units import Quantity, Unit
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .Target import Target

# whether to combine data structures into combined plots
COMBINE_PLOTS = {
    "image": "individual",
    "lightcurve": "combined",
    "spectrum": "individual",
    "sed": "individual",
    "hrd": "combined",
    "powspec": "individual",
    "phasefold": "combined",
    "datatable": "individual",
}

# whether to split plots by target
SPLIT_BY_TARGET = {
    "image": True,
    "lightcurve": True,
    "spectrum": True,
    "sed": True,
    "hrd": False,
    "powspec": True,
    "phasefold": True,
    "datatable": True,
}

SPLIT_BY_SURVEY = {
    "image": True,
    "lightcurve": True,
    "spectrum": True,
    "sed": False,
    "hrd": True,
    "powspec": True,
    "phasefold": True,
    "datatable": False,
}


def manage_inplace(structure: any, inplace: bool):
    if inplace:
        return structure
    else:
        if hasattr(structure, "figure") and structure.figure is not None:
            figure = structure.figure
            structure.figure = None
            struct_copy = copy.deepcopy(structure)
            struct_copy.figure = figure

            return struct_copy
        return copy.deepcopy(structure)


@dataclass
class Container:
    _target_key: str | None = None
    _plot_id: str | None = None

    _data_methods: tuple = ()
    _group_data_methods: dict = field(default_factory=dict)
    _plot_methods: dict = field(default_factory=dict)
    _group_plot_methods: dict = field(default_factory=dict)

    def show(self, show_types: bool = False, show_all: bool = False, **kwargs) -> Self:
        """show(self, show_types = False, show_all = False)
        Prints structure to stdout in a human-readable format.

        Parameters
        ----------
        show_types : bool, optional
            If True, print data types of structure attributes.

        show_all : bool, optional
            If True, do not truncate printing of large iterables.

        Returns
        -------
        ``self``
        """

        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_types, show_all, **kwargs)

        return self

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

    def to_dataframe(self) -> DataFrame:
        """
        Combines all array-like attributes of a structure into a :class:`~pandas.DataFrame`.

        Units are not preserved.

        Returns
        -------
        :class:`~pandas.DataFrame`
        """

        from ..io.structure_io import struct_to_dataframe

        return struct_to_dataframe(self)

    def to_table(self) -> Table:
        """
        Combines all array-like attributes of a structure into a :class:`~astropy.table.Table`, preserving units.

        Returns
        -------
        :class:`~astropy.table.Table`
        """

        from ..io.structure_io import struct_to_table

        return struct_to_table(self)

    def to_hdu(self) -> BinTableHDU:
        from ..io.structure_io import struct_to_hdu

        return struct_to_hdu(self)

    to_hdu.__doc__ = get_docstring("to_hdu", hdu_type="BinTableHDU")

    @classmethod
    def from_dataframe(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs):
        from ..io.structure_io import struct_from_dataframe

        return struct_from_dataframe(cls, target, data, **kwargs)

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: DataFrame, **kwargs):
        from ..io.structure_io import struct_from_table

        return struct_from_table(cls, target, data, **kwargs)
