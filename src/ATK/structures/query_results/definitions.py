from dataclasses import dataclass

import pandas as pd

from .base_definitions import QueryResult


@dataclass
class VizierStruct(QueryResult):
    data: pd.DataFrame | None = None
