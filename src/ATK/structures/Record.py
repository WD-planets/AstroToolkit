from dataclasses import dataclass
from typing import Self

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU
from astropy.table import Table
from pandas import DataFrame

from ..utilities.docstrings import get_docstring
from .structures_core import Container
from .Target import Target


@dataclass(repr=False)
class Record(Container):
    """
    Container for storing rows of a `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue. This object stores metadata describing the origin of the data along with the associated data table.

    |

    """

    #: Alias to `Vizier <https://vizier.cds.unistra.fr/>`_ ``catalogue``.
    survey: str | None = None
    #: `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue ID.
    catalogue: str | None = None
    #: DOC_OVERRIDE
    correction: str | None = None
    #: DOC_OVERRIDE
    search_pos: SkyCoord | None = None
    #: Returned `Vizier <https://vizier.cds.unistra.fr/>`_ table.
    table: Table | None = None

    _required = ["catalogue"]

    def __repr__(self):
        if self.survey:
            return f"<{self.survey} ({self.catalogue}) Record>"
        else:
            return f"<{self.catalogue} Record>"

    def to_hdu(self) -> BinTableHDU:
        # overwrites the default to_hdu method due to simplicity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)

    to_hdu.__doc__ = get_docstring("to_hdu", hdu_type="BinTableHDU")

    def to_table(self) -> Table:
        """
        Converts structure into a :class:`~astropy.table.Table`.
        """

        # overwrites the default to_table method due to simplicity
        if self.table:
            return self.table
        else:
            return Table()

    def to_dataframe(self) -> DataFrame:
        """
        Converts structure into a :class:`~pandas.DataFrame`.
        """

        # overwrites the default to_dataframe method due to simplicity
        if self.table:
            return self.table.to_pandas()
        else:
            return DataFrame()

    @classmethod
    def from_table(cls, target: Target | int | SkyCoord, data: Table, **kwargs) -> Self:
        from ..io.structure_io import struct_from_table

        ctnr = struct_from_table(cls, target, data, **kwargs)

        if not isinstance(data, Table):
            raise ValueError(f"Invalid 'data' type, expected Table, got {type(data)}.")

        ctnr.table = data

        return ctnr

    from_table.__func__.__doc__ = get_docstring("from_table", obj="Record", args=", ".join(f"``{p}``" for p in _required))

    @classmethod
    def from_dataframe(cls: any, target: Target | int | SkyCoord, data: DataFrame, **kwargs) -> Self:
        from ..io.structure_io import struct_from_dataframe

        ctnr = struct_from_dataframe(cls, target, data, **kwargs)

        if not isinstance(data, DataFrame):
            raise ValueError(f"Invalid 'data' type, expected DataFrame, got {type(data)}.")

        ctnr.table = Table.from_pandas(data)

        return ctnr

    from_dataframe.__func__.__doc__ = get_docstring("from_dataframe", obj="Record", args=", ".join(f"``{p}``" for p in _required))
