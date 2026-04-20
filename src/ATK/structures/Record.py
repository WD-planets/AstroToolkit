from dataclasses import dataclass

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io.fits import BinTableHDU
from astropy.table import Table
from pandas import DataFrame

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
    #: Achieved degree of proper motion correction.
    #:
    #: - ``'full'`` = complete 3-dimensional projection on the sky.
    #: - ``'partial'`` = 2-dimensional plane projection.
    #: - ``'none'`` = no correction.
    correction: str | None = None
    #: Position of search at time of execution (i.e. post-correction).
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
        """
        Converts structure into a :class:`~astropy.io.fits.BinTableHDU`.
        """

        # overwrites the default to_hdu method due to simplicity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)

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
    def from_table(cls, target: Target | int | SkyCoord, data: Table, **kwargs):
        from ..io.structure_io import struct_from_table

        ctnr = struct_from_table(cls, target, data, **kwargs)

        if not isinstance(data, Table):
            raise ValueError(f"Invalid 'data' type, expected Table, got {type(data)}.")

        ctnr.table = data

        return ctnr

    @classmethod
    def from_dataframe(cls: any, target: Target | int | SkyCoord, data: pd.DataFrame, **kwargs):
        from ..io.structure_io import struct_from_dataframe

        ctnr = struct_from_dataframe(cls, target, data, **kwargs)

        if not isinstance(data, DataFrame):
            raise ValueError(f"Invalid 'data' type, expected DataFrame, got {type(data)}.")

        ctnr.table = Table.from_pandas(data)

        return ctnr
