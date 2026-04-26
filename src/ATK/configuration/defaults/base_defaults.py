from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_config.ini"

DEFAULTS = {
    "global_settings": {"unit_format": "symbol"},
    "plot_settings": {
        "cache_time": 3600,
        "size": 500,
        "backend": "canvas",
        "toolbars": True,
        "grids": True,
        "titles": True,
        "font": "Helvetica",
        "font_size": 14,
    },
    "query_settings": {"query_radius": 3, "image_size": 30, "default_scale": "arcsec"},
    "overlay_settings": {"crossmatch_radius": 5, "simbad_radius": 3},
    "datapage_settings": {"grid_size": 200, "font": "Helvetica", "font_size": 12},
}


CONFIG_KEY_DEFS = {
    "global_settings": {
        "unit_format": "Sets the default format when printing an astropy :class:`~astropy.units.Quantity`, from (symbol, text)."
    },
    "plot_settings": {
        "cache_time": "Sets the duration for which figures (i.e. HTML pages) are cached to the ``~/.AstroToolkit/cached_figures`` directory in seconds.",
        "size": "Sets the default grid size of figures.",
        "backend": "Sets the bokeh output backend (see `here <https://docs.bokeh.org/en/latest/docs/user_guide/output.html>`_)",
        "toolbars": "Enables/Disables toolbars in figures.",
        "grids": "Enables/Disables background grid in figures.",
        "titles": "Enables/Disables figure titles.",
        "font": "Sets the font to use in all figure text.",
        "font_size": "Sets the font size to use in all figure text.",
    },
    "query_settings": {
        "query_radius": "Sets the default ``radius`` in :func:`~ATK.Tools.query`.",
        "image_size": "Sets the default ``size`` in image queries with :func:`~ATK.Tools.query`.",
        "default_scale": "Sets the default positional scale for ``radius`` and ``size`` in :func:`~ATK.Tools.query`.",
    },
    "overlay_settings": {
        "crossmatch_radius": "Radius to use when matching Gaia detections to non-Gaia detections when performing proper motion correction in image overlays.",
        "simbad_radius": "Radius to use in SIMBAD searches for image overlays (e.g. when clicking on detection markers).",
    },
    "datapage_settings": {
        "grid_size": "Size of grid used to scale all figures in **datapage** plotting.",
        "font": "Font to apply to all text in **datapage** figures.",
        "font_size": "Font size to apply to all text in **datapage** figures.",
    },
}
