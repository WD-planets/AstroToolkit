import copy
from dataclasses import dataclass, field

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity

from ..utilities.docstrings import get_docstring


def skycoord_equality_check(coords, other):
    same_position = False
    if coords.separation(other) < 1 * u.arcsec:
        same_position = True

    if hasattr(other, "obstime") and other.obstime is not None:
        same_time = abs(coords.obstime - other.obstime) < 1e-9 * u.day

        return same_time and same_position

    return same_position


@dataclass
class Target:
    #: Initial (i.e. uncorrected) coordinates of target source.
    #:
    #: If an ``identifier`` is provided, a :class:`~astropy.coordinates.SkyCoord` is generated in the frame and epoch of the chosen astrometric backend survey. Otherwise the input :class:`~astropy.coordinates.SkyCoord` is used directly.
    initial_coords: SkyCoord
    #: Coordinates of target source, corrected for proper motion. Updated with each stage of astrometric correction.
    coords: SkyCoord
    #: Radius of search in which the Target is being utilised.
    radius: Quantity | None = None

    #: Unique identifier from a supported astrometric backend survey (if provided, otherwise ``None``).
    #:
    #: (e.g. a Gaia Source ID).
    identifier: int | None = None
    #: Survey from which ``identifier`` originates (if one was provided, otherwise ``None``).
    survey: str | None = None
    #: Degree of proper motion correction that the :class:`~ATK.Models.Target` can support.
    #:
    #:     - ``'full'`` = complete 3-dimensional projection on the sky.
    #:     - ``'partial'`` = 2-dimensional plane projection.
    #:     - ``'none'`` = no correction.
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
        if isinstance(other, SkyCoord):
            init_coords_match = skycoord_equality_check(self.initial_coords, other)
        else:
            init_coords_match = skycoord_equality_check(self.initial_coords, other.initial_coords)

            identifier_match = self.identifier == other.identifier
            survey_match = self.survey == other.survey
            correction_match = self.correction == other.correction

            matches = [init_coords_match, identifier_match, survey_match, correction_match]

            return all(matches)

        return init_coords_match

    @property
    def _frame(self):
        return self.coords.frame.name

    @property
    def _epoch(self):
        return self.coords.obstime.fits

    @property
    def _initial_frame(self):
        return self.initial_coords.frame.name

    @property
    def _initial_epoch(self):
        return self.initial_coords.obstime.fits

    def show(self, show_types: bool = False, show_all: bool = False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_types, show_all, **kwargs)

    show.__doc__ = get_docstring("show")

    @classmethod
    def from_id(cls, id: int, survey="gaia"):
        """
        Construct a :class:`~ATK.Models.Target` using a unique identifier from a supported astrometric backend survey.

        The resulting :class:`~ATK.Models.Target` is initialised using the coordinate frame and reference epoch of the selected ``survey``. All astrometric parameters (e.g. position and proper motion) are retrieved from the survey’s `VizieR <https://vizier.cds.unistra.fr/>`_ catalogue.

        Parameters
        ----------
        id: int
            Unique source identifier from a supported astrometric backend survey.

            E.g. a Gaia DR3 Source ID

        survey: {'gaia'}, optional
            Survey to use as an astrometric backend.

            Currently only Gaia DR3 (``gaia``) is supported.

        Returns
        -------
        ``Self``
        """

        if survey == "gaia":
            from ..utilities.coordinates import get_gaia_target

            return get_gaia_target(id)
        else:
            raise NotImplementedError("Other astronometric surveys will be added at a later date.")

    @classmethod
    def from_coord(cls, position: SkyCoord):
        """
        Construct a :class:`~ATK.Models.Target` from a position on the sky.

        The provided :class:`~astropy.coordinates.SkyCoord` is transformed to the ICRS
        frame and used to initialise the target. If no observation epoch is defined
        (``obstime``), an epoch of J2000 is assumed.

        Astrometric correction (i.e. propagation between epochs) is only possible
        if the input coordinate includes proper motion information (and optionally
        distance). If these are not provided, no correction can be applied.

        Parameters
        ----------
        position : :class:`~astropy.coordinates.SkyCoord`
            Input sky coordinate defining the target position. May optionally include
            proper motion and distance information.

        Returns
        -------
        ``Self``
        """

        # if no epoch was set, assume J2000
        if not position.obstime:
            j2000 = Time("2000-01-01T00:00:00.000", format="fits")

            position = SkyCoord(position.data, frame=position.frame, obstime=j2000)

        icrs_pos = position.transform_to("icrs")

        return cls(copy.deepcopy(icrs_pos), copy.deepcopy(icrs_pos), None, None, None, "none")
