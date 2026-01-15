import importlib

from astropy.coordinates import SkyCoord

from ..configuration.base_config import BASE_CONFIG
from ..configuration.epoch_config import EPOCH_CONFIG
from ..structures.definitions import (BaseQueryResult, PlottableQueryResult,
                                      QueryResult, Target)
from ..utilities.coordinates import correct_target, prepare_search
from ..utilities.defaults import RETURNS
from ..utilities.mapping import build_map, get_query_result_map


def _normalise_targeting_input(targeting: any):
    """
    Returns a flattened list of base target inputs (int/SkyCoord/Target)
    """

    if isinstance(targeting, (list, tuple)):
        items = []
        for t in targeting:
            # recursive translation
            items.extend(_normalise_targeting_input(t))
        return items

    if isinstance(targeting, SkyCoord) and not targeting.isscalar:
        return list(targeting)

    return [targeting]


def _make_target(obj, astrometric_backend):
    # already a Target
    if isinstance(obj, Target):
        return obj

    # SkyCoord -> Target
    if isinstance(obj, SkyCoord):
        return Target.from_pos(obj)

    # id -> Target
    if isinstance(obj, int):
        return Target.from_id(obj, astrometric_backend)  # can be NULL or EXCEPTION

    raise TypeError(f"Unsupported target type: {type(obj)}")


def setup_targeting(kind: str, targeting, **arguments) -> list[Target]:
    """
    Normalise user targeting input into a list of Targets
    """

    backend = BASE_CONFIG.get("global_settings", "astrometric_backend")

    # normalise into flat list of things that can be turned into Targets
    try:
        items = _normalise_targeting_input(targeting)
    except Exception:
        # nothing should go be able to go wrong but just in case
        raise ValueError("Unexpected error occured in target creation.")

    # construct targets
    targets = []
    for obj in items:
        target = _make_target(obj, backend)
        if target is RETURNS.EXCEPTION:
            return target
        if target is not RETURNS.NULL:
            targets.append(target)

    if not targets:
        return RETURNS.NULL

    return targets


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


def single_target_query(kind: str, target: Target, structure: BaseQueryResult, **arguments):
    """
    Sets up a query on a single target
    """

    module = importlib.import_module(f"ATK.queries.{kind}")
    query_map = build_map(module, "query", suffix="_query")

    # get specific query function (for given survey if multiple are available)
    query_function = query_map[arguments.get("survey")] if len(query_map) > 1 else list(query_map.values())[0]

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

    for ctnr in structure.data:
        ctnr.correction = target.correction
        ctnr.epoch = target.coords.obstime.fits

    return structure


def general_query(kind: str, targets: list[Target | int | SkyCoord], **arguments):
    """
    Sets up a query on multiple targets
    """

    corrected_targets, structure = prepare_search(targets=targets, query_kind=kind, **arguments)

    for target in corrected_targets:
        structure = single_target_query(kind, target, structure, **arguments)

    return structure
