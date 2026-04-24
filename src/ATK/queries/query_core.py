import importlib
import os
import warnings
from types import FunctionType

from astropy.coordinates import SkyCoord
from astropy.io import fits

from ..configuration.base_config import BASE_CONFIG
from ..structures.DataSet import DataSet
from ..structures.Target import Target
from ..Tools.read import read
from ..utilities.coordinates import correct_target, prepare_search
from ..utilities.defaults import RETURNS
from ..utilities.mapping import build_map
from ..utilities.misc import get_package_version
from .checksum import make_cache_key


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
        return Target.from_coord(obj)

    # id -> Target
    if isinstance(obj, int):
        return Target.from_id(obj, astrometric_backend)  # can be NULL or EXCEPTION

    raise TypeError(f"Unsupported target type: {type(obj)}")


def setup_targeting(targeting: any) -> list[Target]:
    """
    Normalise user targeting input into a list of Targets
    """

    backend = BASE_CONFIG._get("global_settings", "astrometric_backend")

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


def _set_results(structure: DataSet, query_result: any) -> DataSet:
    """
    Sets the .data and .exception attributes of an ATK structure based on what was returned from a query
    """

    # an exception was encountered
    if query_result is RETURNS.EXCEPTION:
        structure.exception = True
        return structure

    # no data was returned
    if query_result is RETURNS.NULL:
        return structure

    # data was returned correctly
    if not isinstance(query_result, list):
        raise Exception(f"Unexpected query_result type '{type(query_result)}', expected list.")

    ctnr_types = set(type(ctnr) for ctnr in query_result)
    if len(ctnr_types) > 1:
        # should be impossible
        raise Exception("Query returned more than one container type.")

    structure.data.extend(query_result)

    structure.kind = list(ctnr_types)[0].__name__

    return structure


def set_container_keys(target: Target, query_result: any) -> any:
    if query_result is RETURNS.EXCEPTION:
        return query_result

    if query_result is RETURNS.NULL:
        return query_result

    for ctnr in query_result:
        # e.g. SED gets a correction array instead, don't want to overwrite this
        if hasattr(ctnr, "correction") and ctnr.correction is None:
            ctnr.correction = target.correction
        ctnr._target_key = target._key

    return query_result


def image_requery(query_function: FunctionType, target: Target, structure: DataSet, query_result: list, **arguments) -> tuple[Target, list]:
    from .image.overlays import get_overlay

    image_time = query_result[0].search_pos.obstime
    corrected_pos = correct_target(target, epoch=image_time)

    query_result = query_function(corrected_pos, **arguments)
    if not query_result or query_result is RETURNS.NULL:
        return query_result
    if query_result is RETURNS.EXCEPTION:
        structure.exception = True
        return query_result

    overlay = get_overlay(target, query_result[0], **arguments)

    if overlay is RETURNS.EXCEPTION:
        query_result[0].overlay = None
        structure.exception = True
    else:
        query_result[0].overlay = overlay

    return target, query_result


def single_target_query(kind: str, target: Target, structure: DataSet, **arguments):
    """
    Sets up a query on a single target
    """

    module = importlib.import_module(f".queries.{kind}", package="ATK")
    query_map = build_map(module, "query", suffix="_query")

    # get specific query function (for given survey if multiple are available)
    if len(query_map) > 1:
        query_function = query_map.get(arguments.get("survey"))
        if query_function is None:
            raise ValueError(f"Invalid {kind} survey '{arguments.get('survey')}'. Valid surveys are: {', '.join(query_map.keys())}.")
    else:
        funcs = list(query_map.values())
        if len(funcs):
            query_function = funcs[0]
        else:
            raise Exception("Unexpected query mapping error.")  # shouldn't happen

    # perform query
    query_result = query_function(target, **arguments)
    if query_result is RETURNS.NULL:
        return structure
    if query_result is RETURNS.EXCEPTION:
        structure.exception = True
        return structure

    # perform second query in image queries (at image-corrected position)
    if kind == "image":
        target, query_result = image_requery(query_function, target, structure, query_result, **arguments)

    # set data and exception attributes
    structure = _set_results(structure, query_result)

    # set any additional container keys on each returned container
    query_result = set_container_keys(target, query_result)

    return structure


def recreate_struct(kind: str, targeting: list[Target], **arguments) -> DataSet:
    structure = read(arguments["path"])

    primary_header = fits.open(arguments["path"])[0].header
    if not primary_header.get("ATK_VER"):
        warnings.warn("Could not determine ATK version from local file.")
    else:
        if primary_header["ATK_VER"] != get_package_version():
            warnings.warn(f"ATK version has changed since file '{arguments['path']}' was generated. Query will be re-run and local file will be overwritten.")
            return None

        if primary_header["ATK_EXCEPTION"]:
            warnings.warn(f"Exception was encountered during data retrieval for file '{arguments['path']}'.Query will be re-run and local file will be overwritten.")
            return None

    current_key = make_cache_key(kind, targeting, arguments)

    # print(f"Recreated key:\n{structure._cache_key}\nCurrent key:\n{current_key}")

    if getattr(structure, "_cache_key", None) != current_key:
        warnings.warn("Detected change in query parameters, query will be re-run and local file will be overwritten.")
        return None

    return structure


def general_query(kind: str, targets: list[Target | int | SkyCoord], **arguments):
    """
    Dispatches a query on one or multiple targets
    """

    corrected_targets, structure = prepare_search(targets=targets, query_kind=kind, **arguments)

    if arguments.get("path") and os.path.exists(arguments["path"]):
        # passing corrected targets not actually necessary currently, but might want to check corrected coords in future
        rec_structure = recreate_struct(kind, corrected_targets, **arguments)

        if rec_structure is not None:
            return rec_structure

    for target in corrected_targets:
        structure = single_target_query(kind, target, structure, **arguments)

    return structure
