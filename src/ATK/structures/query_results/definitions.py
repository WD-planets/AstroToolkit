from dataclasses import dataclass

import pandas as pd

from .base_definitions import QueryResult


@dataclass
class VizierStruct(QueryResult):
    data: pd.DataFrame | None = None

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return super().__str__()
