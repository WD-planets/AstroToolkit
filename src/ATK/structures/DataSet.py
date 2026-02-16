from dataclasses import dataclass, field
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from bokeh.plotting import figure as Figure

from .structures_core import COMBINE_PLOTS, SPLIT_BY_TARGET, Container, manage_inplace
from .Target import Target


@dataclass
class BaseDataSet:
    kind: str | None = None
    survey: str | None = None
    targets: list[Target] | None = field(default_factory=list)
    radius: Quantity | None = None
    exception: bool | None = False

    # maps per-Target key to Target
    _key_map: dict[str, Target] = field(init=False, default_factory=dict)
    # maps per-Target alias to per-Target key
    _alias_map: dict[str, str] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self._build_target_maps()

    def __repr__(self):
        return f"<{self.survey} {self.kind} data>"

    def __str__(self):
        return self.__repr__()

    def _build_target_maps(self):
        self._key_map = {t._key: t for t in self.targets}
        self._alias_map = {}
        for t in self.targets:
            for alias in t._aliases:
                self._alias_map[alias] = t._key

    def show(self, show_all_types=False, **kwargs) -> None:
        from ..io.struct_stdout import pprint_structure

        pprint_structure(self, show_all_types, **kwargs)

        return self

    def save(self, path: str | Path = None) -> Path:
        from ..io.files.writing import write_local

        return write_local(self, path)

    def _fetch_by_key(self, key: str):
        return [d for d in self.data if d._target_key == key]

    def fetch_by_id(self, id: int):
        key = self._alias_map.get(f"id:{id}")
        if key is None:
            return []
        return self._fetch_by_key(key)

    def fetch_by_coord(self, coord: SkyCoord, radius: Quantity | None = 3 * u.arcsec):
        for t in self.targets:
            if coord.separation(t.initial_coords) < radius:
                return self._fetch_by_key(t._key)
        return []

    def fetch_by_target(self, target: Target):
        return self._fetch_by_key(target._key)

    @property
    def _fname(self):
        if self.targets[0].identifier:
            fname = f"{self.targets[0].identifier}"
        else:
            fname = f"{self.position.ra.value:.3f}_{self.position.dec.value:.3f}"
        if self.survey:
            fname += f"_{self.survey}"
        if len(self.targets) > 1:
            fname += "_multi"
        fname += f"_ATK{self.kind}.fits"

        return fname

    def apply(self, method: str, *args, inplace=True, **kwargs):
        from .methods.apply import apply_methods

        struct = manage_inplace(self, inplace)

        if len(set(type(ctnr) for ctnr in struct.data)) > 1:
            raise ValueError("DataSet contains more than one kind of Container.")

        struct = apply_methods(struct, method, *args, **kwargs)

        return struct


@dataclass(repr=False)
class DataSet(BaseDataSet):
    data: list = field(default_factory=list)
    figure: Figure | None = None

    _stored_plot_params: dict = field(default_factory=dict)

    @property
    def _title(self):
        return f"ATK {self.kind.upper()}"

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

    def plot(self, **kwargs: any):
        if not self._ctnr_kind:
            return self

        if self._ctnr_kind not in COMBINE_PLOTS:
            raise ValueError(f"{self._ctnr_kind} containers do not support plotting.")

        from ..io.plot_io import plot_data

        self.figure = plot_data(self, **kwargs)
        self._stored_plot_params = kwargs

        return self

    def open(self, fname: Path | str | None = None, **kwargs: any):
        from ..io.plot_io import open as open_html

        open_html(self, fname=fname, **kwargs)

        return self

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
