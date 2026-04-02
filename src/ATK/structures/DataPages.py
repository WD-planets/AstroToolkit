from dataclasses import dataclass, field
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from bokeh.io import output_file
from bokeh.io import save as bokeh_save
from bokeh.layouts import GridBox

from .Target import Target


@dataclass
class DataPages:
    targets: list[Target] = field(default_factory=list)
    figures: list[GridBox] = field(default_factory=list)

    # maps per-Target key to Target
    _key_map: dict[str, Target] = field(init=False, default_factory=dict)
    # maps per-Target alias to per-Target key
    _alias_map: dict[str, str] = field(init=False, default_factory=dict)
    # maps targets to plots
    _plot_map: dict[str, str] = field(default_factory=dict)

    def __repr__(self):
        return f"<{len(self.figures)} DataPages>"

    def __str__(self):
        return self.__repr__()

    def __post_init__(self):
        self._build_target_maps()

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

    def open(self, fname: Path | str | None = None):
        from ..io.plot_io import open_basic as open_html

        if not fname.endswith(".html"):
            fname = f"{fname}.html"

        for index, plot in enumerate(self.figures):
            if fname:
                out_name = f"{fname[:-5]}_{index + 1}.html"
            else:
                out_name = None

            open_html(plot, fname=out_name, title="ATK DATAPAGE")

    def save(self, fname: Path | str):
        if not fname.endswith(".html"):
            fname = f"{fname}.html"

        for index, plot in enumerate(self.figures):
            out_name = f"{fname[:-5]}_{index + 1}.html"

            output_file(out_name)
            bokeh_save(plot, title="ATK DATAPAGE")

    def _open_by_key(self, key: str, fname: Path | str | None = None):
        from ..io.plot_io import open_basic as open_html

        plot_id = self._plot_map[key]
        for plot in self.figures:
            if plot.id == plot_id:
                open_html(plot, fname=fname, prefix="datapage_", title="ATK DATAPAGE")

    def open_by_id(self, id: int, fname: Path | str | None = None):
        key = self._alias_map.get(f"id:{id}")
        if key is None:
            return
        self._open_by_key(key, fname)

    def open_by_coord(self, coord: SkyCoord, fname: Path | str | None = None, radius: Quantity | None = 3 * u.arcsec):
        for t in self.targets:
            if coord.separation(t.initial_coords) < radius:
                self._open_by_key(t._key, fname)

    def open_by_target(self, target: Target, fname: Path | str | None = None):
        self._open_by_key(target._key, fname)

    def _save_by_key(self, key: str, fname: Path | str):
        plot_id = self._plot_map[key]
        for plot in self.figures:
            if plot.id == plot_id:
                output_file(fname)
                bokeh_save(plot, title="ATK DATAPAGE")

    def save_by_id(self, id: int, fname: Path | str):
        key = self._alias_map.get(f"id:{id}")
        if key is None:
            return
        self._save_by_key(key, fname)

    def save_by_coord(self, coord: SkyCoord, fname: Path | str, radius: Quantity | None = 3 * u.arcsec):
        for t in self.targets:
            if coord.separation(t.initial_coords) < radius:
                self._save_by_key(t._key, fname)

    def save_by_target(self, target: Target, fname: Path | str):
        self._save_by_key(target._key, fname)
