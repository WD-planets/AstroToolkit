from dataclasses import dataclass

import astropy.units as u
import pandas
from astropy.coordinates import SkyCoord
from astropy.io.fits import ImageHDU
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs import WCS

from ..configuration.base_config import BASE_CONFIG
from .structures_core import Container

default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")


@dataclass(repr=False)
class Image(Container):
    survey: str | None = None
    correction: str | None = None
    search_pos: SkyCoord | None = None
    band: str | None = None

    size: Quantity | None = None
    epoch: Time | None = None
    hdu: ImageHDU | None = None
    wcs: WCS | None = None
    overlay: pandas.DataFrame | None = None

    _units = {"size": default_scale}

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def to_hdu(self):
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import image_to_hdu

        return image_to_hdu(self)
