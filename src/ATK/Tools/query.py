import warnings

from astropy.coordinates import SkyCoord

from ..queries.arguments import get_query_arguments
from ..queries.checksum import make_cache_key
from ..queries.query_core import general_query, setup_targeting
from ..structures.DataSet import DataSet
from ..utilities.defaults import RETURNS


def query(kind: str, **arguments) -> DataSet:
    """
    Central query function
    """

    # get necessary parameters from config if not given
    arguments = get_query_arguments(kind, arguments)

    targets = arguments.pop("targets", None)

    if targets is None:
        raise ValueError("No targets specified for query.")

    # get flattened list of targets
    targets = setup_targeting(targets)

    if targets is RETURNS.NULL:
        raise ValueError("Query received no targets.")

    if targets is RETURNS.EXCEPTION:
        raise Exception("Failed to generate requested targets, this is likely due to a Vizier fault.")

    # targets may return exception if Vizier is down
    if any(target is RETURNS.EXCEPTION for target in targets):
        warnings.warn("Failed to generate requested targets, this is likely due to a Vizier fault.")

        structure = DataSet[kind](
            kind=kind, survey=arguments.get("survey", None), targets=None, radius=arguments.get("radius", None), exception=True
        )
        return structure

    # disable proper motion correction
    if arguments.get("disable_correction", False):
        for target in targets:
            initial_coords = target.initial_coords
            target.initial_coords = SkyCoord(
                ra=initial_coords.ra, dec=initial_coords.dec, frame=initial_coords.frame, obstime=initial_coords.obstime
            )
            target.coords = initial_coords
            target.identifier = None
            target.survey = None
            target.correction = "none"

    structure = general_query(kind, targets, **arguments)

    # attach cache key
    structure._cache_key = make_cache_key(kind, targets, arguments)

    if structure.exception:
        return structure

    # save if path provided
    if arguments.get("path"):
        structure.store(arguments["path"])

    return structure
