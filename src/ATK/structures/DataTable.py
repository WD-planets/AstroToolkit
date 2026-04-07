from dataclasses import dataclass

import pandas

from .structures_core import Container


@dataclass(repr=False)
class DataTable(Container):
    table: pandas.DataFrame | None = None

    def __repr__(self):
        return f"<{type(self).__name__}>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)
