from dataclasses import dataclass

import pandas

from .structures_core import Container


@dataclass(repr=False)
class Record(Container):
    survey: str | None = None
    catalogue: str | None = None
    correction: str | None = None
    search_pos: str | None = None
    separation: str | None = None

    data: pandas.DataFrame | None = None

    def __repr__(self):
        if self.survey:
            return f"<{self.survey} ({self.catalogue}) Record>"
        else:
            return f"<{self.catalogue} Record>"

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import simple_to_hdu

        return simple_to_hdu(self)
