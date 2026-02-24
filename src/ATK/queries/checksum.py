import hashlib
import json
from dataclasses import asdict, is_dataclass

from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from ..structures.Target import Target


def _normalise(obj):
    """
    Convert complex objects into deterministic, JSON-serializable form.
    """

    if isinstance(obj, Target):
        return {"initial_coords": obj.initial_coords, "identifier": obj.identifier, "survey": obj.survey}

    # needed to stop infinite recursion
    if isinstance(obj, Quantity):
        return {"value": obj.value.tolist() if hasattr(obj.value, "tolist") else obj.value, "unit": str(obj.unit)}

    # ^
    if isinstance(obj, SkyCoord):
        return {"ra": obj.ra.deg, "dec": obj.dec.deg, "frame": obj.frame.name, "obstime": obj.obstime.isot if obj.obstime else None}

    # dict
    if isinstance(obj, dict):
        return {str(k): _normalise(v) for k, v in sorted(obj.items(), key=lambda item: str(item[0]))}

    # set
    if isinstance(obj, set):
        return sorted(_normalise(v) for v in obj)

    # array-like
    if isinstance(obj, (list, tuple)):
        return [_normalise(v) for v in obj]

    # dataclasses
    if is_dataclass(obj):
        return _normalise(asdict(obj))

    # objects with __dict__
    if hasattr(obj, "__dict__"):
        return _normalise(vars(obj))

    # numpy arrays
    if hasattr(obj, "tolist"):
        return obj.tolist()

    return obj


def make_cache_key(kind, targeting, arguments):
    payload = {"kind": kind, "targets": targeting, "arguments": {k: v for k, v in arguments.items() if k != "path"}}

    # print(payload)

    normalized = _normalise(payload)

    # print(normalized)

    serialized = json.dumps(normalized, sort_keys=True, separators=(",", ":"), default=str)

    key = hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    return key
