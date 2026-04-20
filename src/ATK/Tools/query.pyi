from typing import Literal, overload

from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from ..structures.DataSet import DataSet
from ..structures.Target import Target

@overload
def query(
    kind: Literal["vizier"], targets: int | SkyCoord | Target | list[int | SkyCoord | Target], radius: float | Quantity
) -> DataSet: ...
@overload
def query(kind: Literal["image"], targets: int | SkyCoord | Target | list[int | SkyCoord | Target], size: float | Quantity) -> DataSet: ...
