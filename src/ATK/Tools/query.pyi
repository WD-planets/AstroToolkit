from typing import Literal, overload

from astropy.units import Quantity

@overload
def query(kind: Literal["vizier"], radius: float | Quantity): ...
@overload
def query(kind: Literal["image"], size: float | Quantity): ...
