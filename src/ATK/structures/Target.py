import copy
from dataclasses import dataclass, field

from astropy.coordinates import SkyCoord
from astropy.time import Time


@dataclass
class Target:
    initial_coords: SkyCoord
    coords: SkyCoord

    identifier: int | None = None
    survey: str | None = None
    correction: str = "none"

    _key: str = field(init=False)
    _aliases: set[str] = field(default_factory=set, init=False)

    def __post_init__(self):
        id_key = f"id:{self.identifier}"
        coord_key = f"coord:{self.initial_coords.ra.deg:.8f},{self.initial_coords.dec.deg:.8f}"

        if self.identifier is not None:
            self._key = id_key
            self._aliases.add(id_key)
        else:
            self._key = coord_key
        self._aliases.add(coord_key)

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

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

    @classmethod
    def from_id(cls, id: int, survey="gaia"):
        if survey == "gaia":
            from ..utilities.coordinates import get_gaia_target

            return get_gaia_target(id)
        else:
            raise NotImplementedError("Other astronometric surveys will be added at a later date.")

    @classmethod
    def from_pos(cls, position: SkyCoord):
        # if no epoch was set, assume J2000
        if not position.obstime:
            j2000 = Time("2000-01-01T00:00:00.000", format="fits")

            position = SkyCoord(position.data, frame=position.frame, obstime=j2000)

        icrs_pos = position.transform_to("icrs")

        return cls(copy.deepcopy(icrs_pos), copy.deepcopy(icrs_pos), None, None, "none")
