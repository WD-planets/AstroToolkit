import configparser
import os

from ..PackageInfo import SurveyInfo

surveyInfo = SurveyInfo()

overlay_marker_surveys = surveyInfo.magSurveys
keys_to_not_lower = ["query_lightcurve_atlas_username", "query_lightcurve_atlas_password", "font"]
for survey in overlay_marker_surveys:
    keys_to_not_lower.append(f"{survey}_overlay_mag")


class ConfigStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.config_file = files("AstroToolkit.Configuration").joinpath("ATKConfig.ini")
        if not os.path.isfile(self.config_file):
            print("No ATKConfig.ini found. Generating one with default values...")
            self.set_default_config()
            self.write_config()

    @property
    def supported_keys(self):
        keys = []
        for key, val in vars(self).items():
            if key != "config_file":
                keys.append(key)
        return keys

    def set_default_config(self):
        self.enable_notifications = "True"
        self.query_data_radius = "3"
        self.query_bulkdata_radius = "3"
        self.query_lightcurve_radius = "3"
        self.query_spectrum_radius = "3"
        self.query_sed_radius = "3"
        self.query_reddening_radius = "3"
        self.query_image_size = "30"
        self.query_image_overlays = "gaia"
        self.query_image_band = "g"
        self.query_lightcurve_atlas_username = "None"
        self.query_lightcurve_atlas_password = "None"
        self.unit_size = "500"
        self.search_radius = "3"
        self.datapage_search_button_radius = "3"
        self.overlay_simbad_search_radius = "3"
        self.datapage_datatable_radius = "3"
        self.datapage_grid_size = "250"
        self.datapage_font_size = "12"
        self.output_backend = "canvas"
        self.show_toolbars = "True"
        self.show_grids = "True"
        self.show_titles = "True"
        self.font_size = "14"
        self.font = "Helvetica"
        self.overlay_piggyback_radius = "5"

    def read_config(self):
        Config = configparser.ConfigParser()
        Config.read(self.config_file)

        sections = []
        sections.append(Config["global_settings"])
        sections.append(Config["query_settings"])
        sections.append(Config["image_overlay_settings"])
        sections.append(Config["search_settings"])
        sections.append(Config["datapage_settings"])

        sections_strings = [
            "global_settings",
            "query_settings",
            "image_overlay_settings",
            "search_settings",
            "datapage_settings",
        ]

        self.structured_out = {}

        # read config values, apply some transformations from str -> None/bool
        for section, section_str in zip(sections, sections_strings):
            self.structured_out[section_str] = []
            for key, val in section.items():
                if key not in keys_to_not_lower:
                    val = val.lower()
                if val == "none":
                    val = None
                elif val == "true":
                    val = True
                elif val == "false":
                    val = False
                if key == "font":
                    val = val.title()

                self.structured_out[section_str].append({key: val})
                setattr(self, key, val)

    def write_config(self):
        config = configparser.ConfigParser()

        # take current config values and cast all to lowercase strings
        for key, val in vars(self).items():
            if key not in ["config_file", "structured_out"]:
                val = str(val)
                if key not in keys_to_not_lower:
                    val = val.lower()
                setattr(self, key, val)

        config.add_section("global_settings")
        config.set("global_settings", "enable_notifications", self.enable_notifications)
        config.set("global_settings", "unit_size", self.unit_size)
        config.set("global_settings", "output_backend", self.output_backend)
        config.set("global_settings", "show_toolbars", self.show_toolbars)
        config.set("global_settings", "show_grids", self.show_grids)
        config.set("global_settings", "show_titles", self.show_titles)
        config.set("global_settings", "font_size", self.font_size)
        config.set("global_settings", "font", self.font)

        config.add_section("query_settings")
        config.set("query_settings", "query_data_radius", self.query_data_radius)
        config.set("query_settings", "query_bulkdata_radius", self.query_bulkdata_radius)
        config.set("query_settings", "query_lightcurve_radius", self.query_lightcurve_radius)
        config.set("query_settings", "query_spectrum_radius", self.query_spectrum_radius)
        config.set("query_settings", "query_sed_radius", self.query_sed_radius)
        config.set("query_settings", "query_reddening_radius", self.query_reddening_radius)
        config.set("query_settings", "query_image_size", self.query_image_size)
        config.set("query_settings", "query_image_overlays", self.query_image_overlays)
        config.set("query_settings", "query_image_band", self.query_image_band)
        config.set("query_settings", "query_lightcurve_atlas_username", self.query_lightcurve_atlas_username)
        config.set("query_settings", "query_lightcurve_atlas_password", self.query_lightcurve_atlas_password)

        config.add_section("image_overlay_settings")
        config.set("image_overlay_settings", "overlay_piggyback_radius", self.overlay_piggyback_radius)
        config.set("image_overlay_settings", "overlay_simbad_search_radius", self.overlay_simbad_search_radius)

        config.add_section("search_settings")
        config.set("search_settings", "search_radius", self.search_radius)

        config.add_section("datapage_settings")
        config.set("datapage_settings", "datapage_search_button_radius", self.datapage_search_button_radius)
        config.set("datapage_settings", "datapage_datatable_radius", self.datapage_datatable_radius)
        config.set("datapage_settings", "datapage_font_size", self.datapage_font_size)
        config.set("datapage_settings", "datapage_grid_size", self.datapage_grid_size)

        with open(self.config_file, "w") as file:
            config.write(file)

    def edit_config(self, key, value):
        self.read_config()
        if not isinstance(value, str):
            try:
                value = str(value)
            except:
                raise Exception("ATKConfig.ini options must be strings.")
        setattr(self, key, value)
        self.write_config()

        self.output_config()

    def output_config(self):
        self.read_config()
        for label, section in self.structured_out.items():
            print(f"[{label}]")
            for entry in section:
                for key, val in entry.items():
                    print(f"{key} = {val}")
            print()
