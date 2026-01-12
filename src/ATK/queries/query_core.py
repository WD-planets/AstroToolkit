import copy
import importlib

from astropy.coordinates import SkyCoord

from ..configuration.base_config import BASE_CONFIG
from ..configuration.epoch_config import EPOCH_CONFIG
from ..structures.definitions import PlottableQueryResult, QueryResult, Target
from ..utilities.coordinates import correct_target, prepare_search
from ..utilities.defaults import RETURNS
from ..utilities.mapping import build_map, get_query_result_map

MULTI_QUERY_KEEP_ATTRS = ["kind", "survey", "positions", "identifiers", "epoch", "frame", "exception", "data", "figure"]


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


def single_target_query(kind: str, target: Target, **arguments):
    """
    Sets up a query on a single target
    """

    # if Vizier is down, a target may fail to generate from an ID
    if target is RETURNS.EXCEPTION:
        structure_map = get_query_result_map()
        structure = structure_map[kind](
            kind=kind,
            survey=arguments.get("survey", None),
            position=None,
            identifier=None,
            radius=arguments.get("radius", None),
            frame=None,
            epoch=None,
            correction="none",
            exception=True,
        )
        return structure

    module = importlib.import_module(f"ATK.queries.{kind}")
    query_map = build_map(module, "query", suffix="_query")

    # get specific query function (for given survey if multiple are available)
    query_function = query_map[arguments.get("survey")] if len(query_map) > 1 else list(query_map.values())[0]

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

    survey = arguments.get("survey", None)
    if survey and survey not in EPOCH_CONFIG.get_section_by_query_kind(kind):
        structure.correction = "none"
    else:
        structure.correction = target.correction

    return structure


def multiple_target_query(kind: str, targets: list[Target | int | SkyCoord], **arguments):
    """
    Sets up a query on multiple targets
    """

    structures = []
    for target in targets:
        structure = single_target_query(kind, target, **arguments)
        structures.append(structure)

    final_structure = copy.deepcopy(structures[0])

    final_structure.data = []
    final_structure.identifiers = "multiple"
    final_structure.positions = "multiple"

    for attr in vars(final_structure):
        if attr not in MULTI_QUERY_KEEP_ATTRS:
            setattr(final_structure, attr, None)

    for structure in structures:
        final_structure.data += structure.data

    return final_structure
