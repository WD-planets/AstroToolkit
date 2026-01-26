import os
import warnings

from astropy.coordinates import SkyCoord

from ..queries.arguments import get_query_arguments
from ..queries.query_core import general_query, setup_targeting
from ..structures.DataSet import DataSet
from ..utilities.defaults import RETURNS
from .read import read


def query(kind: str, **arguments) -> DataSet:
    """
    Central query function
    """

    # get necessary parameters from config if not given
    arguments = get_query_arguments(kind, arguments)

    if arguments.get("path") and os.path.exists(arguments["path"]):
        structure = read(arguments["path"])
        return structure

    target, targets = arguments.pop("target", None), arguments.pop("targets", None)

    if target is None == targets is None:
        raise ValueError("Exactly one of 'target' and 'targets' is required for all queries.")

    # get targets that were provided
    targeting = target if target is not None else targets

    # get flattened list of targets
    targets = setup_targeting(kind, targeting)

    # targets may return exception if Vizier is down
    if any(target is RETURNS.EXCEPTION for target in targets):
        warnings.warn("Failed to generate requested targets, this is likely due to a Vizier fault.")

        structure = DataSet[kind](kind=kind, survey=arguments.get("survey", None), targets=None, radius=arguments.get("radius", None), exception=True)
        return structure

    # disable proper motion correction
    if arguments.get("disable_correction", False):
        for target in targets:
            initial_coords = target.initial_coords
            target.initial_coords = SkyCoord(ra=initial_coords.ra, dec=initial_coords.dec, frame=initial_coords.frame, obstime=initial_coords.obstime)
            target.coords = initial_coords
            target.identifier = None
            target.survey = None
            target.correction = "none"

    return general_query(kind, targets, **arguments)
