from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ....structures.Lightcurve import Lightcurve

from dataclasses import dataclass, field
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from bokeh.plotting import figure as Figure

from .structures_core import COMBINE_PLOTS, SPLIT_BY_SURVEY, SPLIT_BY_TARGET, Container, manage_inplace
from .Target import Target


@dataclass
class DataSet:
    kind: str | None = None
    targets: list[Target] | None = field(default_factory=list)
    exception: bool | None = False

    data: list = field(default_factory=list)
    figure: Figure | None = None

    _stored_plot_params: dict = field(default_factory=dict)
    _plotted_keys: list = field(default_factory=list)

    # maps per-Target key to Target
    _key_map: dict[str, Target] = field(init=False, default_factory=dict)
    # maps per-Target alias to per-Target key
    _alias_map: dict[str, str] = field(init=False, default_factory=dict)

    _cache_key: str | None = None

    def __post_init__(self):
        self._build_target_maps()

    def __repr__(self):
        return f"<{self.kind} DataSet>"

    def __str__(self):
        return self.__repr__()

    def _build_target_maps(self):
        self._key_map = {t._key: t for t in self.targets}
        self._alias_map = {}
        for t in self.targets:
            for alias in t._aliases:
                self._alias_map[alias] = t._key

    def show(self, show_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_types, **kwargs)

        return self

    def store(self, path: str | Path = None) -> Path:
        from ..io.files.writing import write_local

        write_local(self, path)

        return self

    def plot(self, **kwargs: any):
        if not self._ctnr_kind:
            return self

        if self._ctnr_kind not in COMBINE_PLOTS:
            raise ValueError(f"{self._ctnr_kind} containers do not support plotting.")

        from ..io.plot_io import plot_data

        keys = []
        for target in self.targets:
            keys.append(target._key)

        self.figure = plot_data(self, **kwargs)
        self._stored_plot_params = kwargs
        self._plotted_keys = keys

        return self

    def open(self, path: Path | str | None = None, **kwargs: any):
        from ..io.plot_io import open as open_html

        keys = []
        for target in self.targets:
            keys.append(target._key)

        open_html(self, fname=path, keys=keys, **kwargs)

        return self

    def save(self, path: Path | str, **kwargs: any):
        from ..io.plot_io import save

        keys = []
        for target in self.targets:
            keys.append(target._key)

        save(self, fname=path, keys=keys, **kwargs)

        return self

    def apply(self, method: str, *args, inplace=True, **kwargs):
        from .methods.apply import apply_methods

        struct = manage_inplace(self, inplace)

        if len(set(type(ctnr) for ctnr in struct.data)) > 1:
            raise ValueError("DataSet contains more than one kind of Container.")

        struct = apply_methods(struct, method, *args, **kwargs)

        return struct

    # ====================
    # Multi-Target Methods
    # ====================

    def _fetch_by_key(self, key: str):
        return [ctnr for ctnr in self.data if ctnr._target_key == key].copy()

    def split(self, targets: any, inplace: bool = True) -> DataSet:
        from ..queries.query_core import setup_targeting

        targets = setup_targeting(targets)
        struct = manage_inplace(self, inplace)

        ctnrs, out_targets = [], []
        for t in targets:
            ctnrs.extend(struct._fetch_by_key(t._key))

        for t in targets:
            for s_t in struct.targets:
                if t._key == s_t._key:
                    out_targets.append(s_t)

        struct.data = ctnrs
        struct.targets = out_targets

        struct._build_target_maps()

        return struct

    def merge(self, dataset: DataSet, inplace: bool = True) -> DataSet:
        struct = manage_inplace(self, inplace)

        struct.data = struct.data + dataset.data
        for target in dataset.targets:
            struct_keys = [t._key for t in struct.targets]
            if target._key not in struct_keys:
                struct.targets.append(target)

        struct._build_target_maps()

        return struct

    # Other Stuff
    # ===========

    @property
    def _title(self):
        return f"ATK {self._ctnr_kind.upper()}"

    @property
    def _ctnr_kind(self):
        ctnr_kinds = [type(ctnr).__name__.lower() for ctnr in self.data]
        if len(set(ctnr_kinds)) > 1:
            raise ValueError("DataSet contains multiple container types.")
        if not ctnr_kinds:
            return

        return ctnr_kinds[0]

    @property
    def _plot_method(self):
        if self._ctnr_kind not in COMBINE_PLOTS:
            return

        return COMBINE_PLOTS[self._ctnr_kind]

    @property
    def _split_by_target(self):
        if self._ctnr_kind not in SPLIT_BY_TARGET:
            return

        return SPLIT_BY_TARGET[self._ctnr_kind]

    @property
    def _split_by_survey(self):
        if self._ctnr_kind not in SPLIT_BY_SURVEY:
            return

        return SPLIT_BY_SURVEY[self._ctnr_kind]

    @classmethod
    def from_target(cls, kind: str, target: Target | int | SkyCoord, radius: float | Quantity | None = None, survey: str = None):
        from ..queries.query_core import setup_targeting

        # kind, targets, survey, radius, exception, data, figure

        targets = setup_targeting(target)

        return cls(kind=kind, targets=targets, survey=survey, radius=radius, exception=False)

    def add(self, data: Container):
        if self.data and self._ctnr_kind != type(data).__name__.lower():
            raise ValueError(f"Cannot add container of type '{type(data).__name__.lower()}' to DataSet containing {self._ctnr_kind} data.")

        self.data.append(data)
