from dataclasses import dataclass

from astropy.coordinates import SkyCoord
from astropy.time import Time


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
