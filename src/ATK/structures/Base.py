import copy
import types
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
    """
    Base class for all data containers.

    See Also
    --------
    :class:`~ATK.Models.Record`
    :class:`~ATK.Models.Image`
    :class:`~ATK.Models.Lightcurve`
    :class:`~ATK.Models.Powspec`
    :class:`~ATK.Models.Spectrum`
    :class:`~ATK.Models.SED`
    :class:`~ATK.Models.HRD`
    :class:`~ATK.Models.DataTable`
    """

    _target_key: str | None = None
    _plot_id: str | None = None

    def show(self, show_types: bool = False, show_all: bool = False, **kwargs) -> Self:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_types, show_all, **kwargs)

        return self

    show.__doc__ = get_docstring("show")

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


# IO Stuff
# --------


def _clone_classmethod(method):
    func = method.__func__
    return types.FunctionType(func.__code__, func.__globals__, name=func.__name__, argdefs=func.__defaults__, closure=func.__closure__)


def _clone_function(func):
    return types.FunctionType(func.__code__, func.__globals__, name=func.__name__, argdefs=func.__defaults__, closure=func.__closure__)


class DataFrameIOMixin:
    """
    Base class for :class:`~pandas.DataFrame` operations.

    See Also
    --------
    :class:`~ATK.Models.Lightcurve`
    :class:`~ATK.Models.Powspec`
    :class:`~ATK.Models.Spectrum`
    :class:`~ATK.Models.SED`
    :class:`~ATK.Models.HRD`
    """

    @classmethod
    def from_dataframe(cls, target: Target | SkyCoord | int, data: DataFrame, **kwargs):
        """
        Handles :class:`~pandas.DataFrame` to container IO.
        """

        from ..io.structure_io import struct_from_dataframe

        return struct_from_dataframe(cls, target, data, **kwargs)

    def to_dataframe(self) -> DataFrame:
        """
        Handles container to :class:`~pandas.DataFrame` IO.
        """

        from ..io.structure_io import struct_to_dataframe

        return struct_to_dataframe(self)

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        new_func = _clone_classmethod(DataFrameIOMixin.from_dataframe)
        new_func.__doc__ = get_docstring("from_dataframe", obj=cls.__name__, args=", ".join(f"``{p}``" for p in cls._required))
        cls.from_dataframe = classmethod(new_func)

        new_func = _clone_function(DataFrameIOMixin.to_dataframe)
        new_func.__doc__ = get_docstring("to_dataframe")
        cls.to_table = new_func


class TableIOMixin:
    """
    Base class for :class:`~astropy.table.Table` operations.

    See Also
    --------
    :class:`~ATK.Models.Lightcurve`
    :class:`~ATK.Models.Powspec`
    :class:`~ATK.Models.Spectrum`
    :class:`~ATK.Models.SED`
    :class:`~ATK.Models.HRD`
    """

    @classmethod
    def from_table(cls, target: Target | SkyCoord | int, data: DataFrame, **kwargs):
        """
        Handles :class:`~astropy.table.Table` to container IO.
        """

        from ..io.structure_io import struct_from_table

        return struct_from_table(cls, target, data, **kwargs)

    def to_table(self) -> Table:
        """
        Handles container to :class:`~astropy.table.Table` IO.
        """

        from ..io.structure_io import struct_to_table

        return struct_to_table(self)

    def __init_subclass__(cls, **kwargs):
        """
        Attach dynamic docstring to each subclass
        """

        super().__init_subclass__(**kwargs)

        # need to create clone of class method so docstring doesn't get overriden with each subclass
        new_func = _clone_classmethod(TableIOMixin.from_table)
        new_func.__doc__ = get_docstring("from_table", obj=cls.__name__, args=", ".join(f"``{p}``" for p in cls._required))
        cls.from_table = classmethod(new_func)

        new_func = _clone_function(TableIOMixin.to_table)
        new_func.__doc__ = get_docstring("to_table")
        cls.to_table = new_func


class FITSIOMixin:
    """
    Base class for FITS operations.

    See Also
    --------
    :class:`~ATK.Models.Image`
    :class:`~ATK.Models.Lightcurve`
    :class:`~ATK.Models.Powspec`
    :class:`~ATK.Models.Spectrum`
    :class:`~ATK.Models.SED`
    :class:`~ATK.Models.HRD`
    """

    def to_hdu(self) -> BinTableHDU:
        """
        Handles container to :class:`~astropy.io.fits.BinTableHDU` IO.
        """

        from ..io.structure_io import struct_to_hdu

        return struct_to_hdu(self)

    to_hdu.__doc__ = get_docstring("to_hdu", hdu_type="BinTableHDU")
