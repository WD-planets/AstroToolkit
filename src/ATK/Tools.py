from pathlib import Path

from astropy.coordinates import SkyCoord

from .io.files.read import read_local
from .queries.arguments import get_query_arguments
from .queries.query_core import (multiple_target_query, setup_targeting,
                                 single_target_query)
from .structures.definitions import PlottableQueryResult, QueryResult


def query(kind: str, **arguments) -> QueryResult | PlottableQueryResult:
    """
    Central query function
    """

    # get necessary parameters from config if not given
    arguments = get_query_arguments(kind, arguments)

    target, targets = arguments.pop("target", None), arguments.pop("targets", None)

    if target is None == targets is None:
        raise ValueError("Exactly one of 'target' and 'targets' is required for all queries.")

    # get targets that were provided
    targeting = target if target is not None else targets

    targets = setup_targeting(kind, targeting)

    # disable proper motion correction
    if arguments.get("disable_corrections", False):
        for target in targets:
            coords = target.coords
            target.coords = SkyCoord(ra=coords.ra, dec=coords.dec, frame=coords.frame, obstime=coords.obstime)
            target.identifier = None
            target.survey = None
            target.correction = "none"

    if len(targets) > 1:
        return multiple_target_query(kind, targets, **arguments)
    else:
        return single_target_query(kind, targets[0], **arguments)


def read(path: str | Path):
    """
    Reads a local ATK fits file to recreate the original data structure
    """

    return read_local(path)
