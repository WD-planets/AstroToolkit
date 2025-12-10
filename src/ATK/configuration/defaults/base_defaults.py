from pathlib import Path

PATH = Path.home() / ".AstroToolkit" / "ATK_config.ini"

DEFAULTS = {
    "global_settings": {"notifications": True},
    "plot_settings": {
        "cache_time": 3600,
        "size": 500,
        "backend": "canvas",
        "toolbars": True,
        "grids": True,
        "titles": True,
        "font_size": 14,
        "font": "Helvetica",
    },
    "query_settings": {"query_radius": 3, "image_size": 30, "image_overlays": "gaia", "image_band": "g"},
    "overlay_settings": {"piggyback_radius": 5, "search_radius": 3},
    "search_settings": {"search_radius": 3},
    "datapage_settings": {"search_button_radius": 3, "datatable_radius": 3, "font_size": 12, "grid_size": 250},
}
