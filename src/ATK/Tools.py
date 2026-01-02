import importlib
from pathlib import Path

from astropy.coordinates import SkyCoord

from .configuration.epoch_config import EPOCH_CONFIG
from .io.files.read import read_local
from .queries.arguments import get_query_arguments
from .structures.definitions import PlottableQueryResult, QueryResult, Target
from .utilities.coordinates import check_target, correct_target, prepare_search
from .utilities.defaults import RETURNS
from .utilities.mapping import build_map, get_query_result_map


def _set_results(structure: QueryResult | PlottableQueryResult, query_result: any) -> QueryResult | PlottableQueryResult:
    """
    Sets the .data and .exception attributes of an ATK structure based on what was returned from a query
    """

    # an exception was encountered
    if query_result is RETURNS.EXCEPTION:
        structure.data = []
        structure.exception = True
        return structure

    # no data was returned
    if query_result is RETURNS.NULL:
        structure.data = []
        return structure

    # data was returned correctly
    if isinstance(query_result, list):
        structure.data += query_result
    else:
        structure.data.append(query_result)

    return structure


def query(kind: str, target: Target | SkyCoord | int, **kwargs) -> QueryResult | PlottableQueryResult:
    """
    Central query function
    """

    target = check_target(target)

    # if no target was made (e.g. because Vizier is down when using Target.from_id)
    if target is RETURNS.EXCEPTION:
        structure_map = get_query_result_map()
        structure = structure_map[kind](
            kind=kind,
            survey=kwargs.get("survey", None),
            position=None,
            identifier=None,
            radius=kwargs.get("radius", None),
            frame=None,
            epoch=None,
            correction="none",
            exception=True,
        )

        return structure

    # disables any proper motion correction by creating a simplified target with no proper motion information
    if kwargs.get("disable_corrections", False):
        coords = target.coords
        reduced_coords = SkyCoord(ra=coords.ra, dec=coords.dec, frame=coords.frame, obstime=coords.obstime)
        target.coords = reduced_coords
        target.identifier = None
        target.survey = None
        target.correction = "none"

    module = importlib.import_module(f"ATK.queries.{kind}")
    query_map = build_map(module, "query", suffix="_query")

    # get necessary parameters from config if not given
    arguments = get_query_arguments(kind, kwargs)

    # get specific query function (for given survey if multiple are available)
    query_function = query_map[kwargs.get("survey")] if len(query_map) > 1 else list(query_map.values())[0]

    # get search position and structure
    search_pos, structure = prepare_search(target=target, query_kind=kind, **arguments)
    if not search_pos:
        return structure

    # perform query
    query_result = query_function(target, **arguments)

    # set data and exception attributes
    structure = _set_results(structure, query_result)

    # perform second query in image queries (at image-corrected position)
    if kind == "image" and structure.data:
        from .queries.image.overlays import get_overlay

        image_time = query_result[0].focus.obstime
        corrected_pos = correct_target(search_pos, epoch=image_time)

        query_result = query_function(corrected_pos, **arguments)

        # delete initial image
        structure.data = []
        # structure.epoch = image_time
        structure = _set_results(structure, query_result)

        overlay = get_overlay(target, structure.data[0], **arguments)

        if overlay is RETURNS.EXCEPTION:
            structure.data[0].overlay = None
            structure.exception = True
        else:
            structure.data[0].overlay = overlay

    survey = kwargs.get("survey", None)
    if survey and survey not in EPOCH_CONFIG.get_section_by_query_kind(kind):
        structure.correction = "none"
    else:
        structure.correction = target.correction

    return structure


def read(path: str | Path):
    """
    Reads a local ATK fits file to recreate the original data structure
    """

    return read_local(path)
