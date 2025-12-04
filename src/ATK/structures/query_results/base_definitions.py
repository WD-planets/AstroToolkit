from dataclasses import dataclass
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.time import Time

from ...io.files.writing import write_structure
from ...io.struct_stdout import pprint_structure


@dataclass
class QueryResult:
    kind: str | None = None
    survey: str | None = None
    radius: float | None = None
    source: int | None = None
    position: SkyCoord | None = None
    epoch: Time | None = None
    frame: str | None = None
    correction: str | None = None
    exception: bool | None = False

    def show(self, show_all_types=False):
        pprint_structure(self, show_all_types)

    def save(self, path: str | Path = None):
        write_structure(self, path)

    def __repr__(self):
        return f"{self.survey} {self.kind} data"

    def __str__(self):
        return self.__repr__()
