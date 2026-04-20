import copy
from dataclasses import dataclass, field

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity


def skycoord_equality_check(self, other):
    same_position = (
        self.frame.name == other.frame.name and np.isclose(self.ra.deg, other.ra.deg) and np.isclose(self.dec.deg, other.dec.deg)
    )
    same_time = abs(self.obstime - other.obstime) < 1e-9 * u.day

    return same_time and same_position


@dataclass
class Target:
    initial_coords: SkyCoord
    coords: SkyCoord
    radius: Quantity | None = None

    identifier: int | None = None
    survey: str | None = None
    correction: str = "none"

    _key: str = field(init=False)
    _aliases: set[str] = field(default_factory=set, init=False)

    def __repr__(self):
        from ..io.struct_stdout import format_target

        return f"<{format_target(self)} {type(self).__name__}>"

    def __post_init__(self):
        id_key = f"id:{self.identifier}"
        coord_key = f"coord:{self.initial_coords.ra.deg:.8f},{self.initial_coords.dec.deg:.8f}"

        if self.identifier is not None:
            self._key = id_key
            self._aliases.add(id_key)
        else:
            self._key = coord_key
        self._aliases.add(coord_key)

    def __eq__(self, other):
        init_coords_match = skycoord_equality_check(self.initial_coords, other.initial_coords)
        identifier_match = self.identifier == other.identifier
        survey_match = self.survey == other.survey
        correction_match = self.correction == other.correction

        matches = [init_coords_match, identifier_match, survey_match, correction_match]

        return all(matches)

    @property
    def frame(self):
        return self.coords.frame.name

    @property
    def epoch(self):
        return self.coords.obstime.fits

    @property
    def initial_frame(self):
        return self.initial_coords.frame.name

    @property
    def initial_epoch(self):
        return self.initial_coords.obstime.fits

    def show(self, show_types: bool = False, show_all: bool = False, **kwargs) -> None:
        """show(self, show_types = False, show_all = False)
        Prints structure to stdout in a human-readable format.

        Parameters
        ----------
        show_types : bool, optional
            If True, print data types of structure attributes.

            Default is ``False``

        show_all : bool, optional
            If True, do not truncate printing of large iterables.

            Default is ``False``.
        """

        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_types, show_all, **kwargs)

    @classmethod
    def from_id(cls, id: int, survey="gaia"):
        if survey == "gaia":
            from ..utilities.coordinates import get_gaia_target

            return get_gaia_target(id)
        else:
            raise NotImplementedError("Other astronometric surveys will be added at a later date.")

    @classmethod
    def from_coord(cls, position: SkyCoord):
        # if no epoch was set, assume J2000
        if not position.obstime:
            j2000 = Time("2000-01-01T00:00:00.000", format="fits")

            position = SkyCoord(position.data, frame=position.frame, obstime=j2000)

        icrs_pos = position.transform_to("icrs")

        return cls(copy.deepcopy(icrs_pos), copy.deepcopy(icrs_pos), None, None, None, "none")
