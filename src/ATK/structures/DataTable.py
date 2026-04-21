from dataclasses import dataclass
from typing import Self

from astropy.coordinates import SkyCoord
from astropy.table import Table
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .structures_core import Container
from .Target import Target


@dataclass(repr=False)
class DataTable(Container):
    """
    Container for storing target-matched data from a subset of one or multiple `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue.
    """

    #: DOC_OVERRIDE
    table: Table | None = None

    def __repr__(self):
        return f"<{type(self).__name__}>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)

    to_hdu.__doc__ = get_docstring("to_hdu", hdu_type="BinTableHDU")

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: Table, **kwargs) -> Self:
        from ..io.structure_io import struct_from_table

        ctnr = struct_from_table(cls, target, data, **kwargs)

        if not isinstance(data, Table):
            raise ValueError(f"Invalid 'data' type, expected Table, got {type(data)}.")

        ctnr.table = data

        return ctnr

    from_table.__func__.__doc__ = get_docstring("from_table", obj="DataTable", args="None")

    @classmethod
    def from_dataframe(cls: any, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        from ..io.structure_io import struct_from_dataframe

        ctnr = struct_from_dataframe(cls, target, data, **kwargs)

        if not isinstance(data, DataFrame):
            raise ValueError(f"Invalid 'data' type, expected DataFrame, got {type(data)}.")

        ctnr.table = Table.from_pandas(data)

        return ctnr

    from_dataframe.__func__.__doc__ = get_docstring("from_dataframe", obj="DataTable", args="None")
