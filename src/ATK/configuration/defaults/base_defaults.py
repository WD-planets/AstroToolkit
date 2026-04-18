from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_config.ini"

DEFAULTS = {
    "global_settings": {"astrometric_backend": "gaia", "notifications": True, "unit_format": "symbol"},
    "plot_settings": {
        "cache_time": 3600,
        "size": 500,
        "backend": "canvas",
        "toolbars": True,
        "grids": True,
        "titles": False,
        "font": "Helvetica",
        "font_size": 14,
    },
    "query_settings": {"query_radius": 3, "image_size": 30, "default_scale": "arcsec"},
    "overlay_settings": {"piggyback_radius": 5, "simbad_radius": 3},
    "datapage_settings": {"grid_size": 200, "font": "Helvetica", "font_size": 12},
}
