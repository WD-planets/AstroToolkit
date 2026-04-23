from dataclasses import dataclass

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import ImageHDU
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs import WCS
from pandas import DataFrame

from ..configuration.base_config import BASE_CONFIG
from ..utilities.docstrings import get_docstring
from .structures_core import Container

default_scale = BASE_CONFIG._get("query_settings", "default_scale")
try:
    default_unit = u.Unit(default_scale)
except ValueError:
    raise Exception(f"Invalid default_unit in config '{default_scale}'.")

PLOT_PARAMS = {
    "cmap": [
        "{'viridis', 'grey', 'false_colour'}, optional",
        """
Sets the image colour map.

``false_colour`` maps the wavelength of the image filter into a single RGB colour.

Default is ``'viridis'``.
""",
    ],
    "relative_axes": [
        "bool, optional",
        """
Sets the format of the coordinate axes.

``True``: coordinate axes display distance from the image centre.

``False``: coordinate axes display positions on the sky.

Default is ``True``.
""",
    ],
}


@dataclass(repr=False)
class Image(Container):
    """
    Container for storing image data. This object stores both data and relevant metadata.
    """

    #: DOC_OVERRIDE
    survey: str | None = None
    #: DOC_OVERRIDE
    correction: str | None = None
    #: DOC_OVERRIDE
    search_pos: SkyCoord | None = None
    #: DOC_OVERRIDE
    band: str | None = None

    #: Size of stored image.
    size: Quantity | None = None
    #: Epoch of stored image.
    epoch: Time | None = None
    #: Raw image data and header.
    hdu: ImageHDU | None = None
    #: World coordinate system of stored image.
    wcs: WCS | None = None
    #: DOC_OVERRIDE
    overlay: DataFrame | None = None

    _units = {"size": default_scale}

    _plot_params = PLOT_PARAMS

    def __repr__(self):
        return f"<{self.survey} {self.band}-band {type(self).__name__}>"

    def __post_init__(self):
        for attr, unit in self._units.items():
            val = getattr(self, attr, None)
            if val is None:
                continue

            if not isinstance(val, Quantity):
                setattr(self, attr, val * unit)

    def to_hdu(self) -> ImageHDU:
        # overwrites the default to_hdu method due to complexity
        from ..io.structure_io import image_to_hdu

        return image_to_hdu(self)

    to_hdu.__doc__ = get_docstring("to_hdu", hdu_type="ImageHDU")
