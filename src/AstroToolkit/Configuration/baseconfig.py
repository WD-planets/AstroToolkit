import configparser
import os


class ConfigStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.config_file = files("AstroToolkit.Configuration").joinpath("ATKConfig.ini")
        if not os.path.isfile(self.config_file):
            print("No ATKConfig.ini found. Making new one with default values\n")
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
        from ..Data.dataquery import SurveyInfo

        overlay_info = SurveyInfo().overlay_param_names

        self.enable_notifications = "True"
        self.query_data_radius = "3"
        self.query_phot_radius = "3"
        self.query_bulkphot_radius = "3"
        self.query_lightcurve_radius = "3"
        self.query_spectrum_radius = "3"
        self.query_sed_radius = "3"
        self.query_reddening_radius = "5"
        self.query_image_size = "30"
        self.query_image_overlays = "gaia"
        self.query_image_band = "g"
        self.query_lightcurve_atlas_username = "None"
        self.query_lightcurve_atlas_password = "None"
        self.unit_size = "400"
        self.search_radius = "3"
        self.datapage_search_button_radius = "3"
        self.gaia_overlay_mag = overlay_info["gaia"]["default_overlay_mag"]
        self.galex_overlay_mag = overlay_info["galex"]["default_overlay_mag"]
        self.wise_overlay_mag = overlay_info["wise"]["default_overlay_mag"]
        self.sdss_overlay_mag = overlay_info["sdss"]["default_overlay_mag"]
        self.twomass_overlay_mag = overlay_info["twomass"]["default_overlay_mag"]
        self.skymapper_overlay_mag = overlay_info["skymapper"]["default_overlay_mag"]
        self.panstarrs_overlay_mag = overlay_info["panstarrs"]["default_overlay_mag"]
        self.datapage_datatable_radius = "3"
        self.datapage_grid_size = "250"
        self.output_backend = "canvas"
        self.font_size = "12"
        self.font = "Times New Roman"
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

        for section in sections:
            for key, val in section.items():
                setattr(self, key, val)

    def write_config(self):
        config = configparser.ConfigParser()

        config.add_section("global_settings")
        config.set("global_settings", "enable_notifications", self.enable_notifications)
        config.set("global_settings", "unit_size", self.unit_size)
        config.set("global_settings", "output_backend", self.output_backend)
        config.set("global_settings", "font_size", self.font_size)
        config.set("global_settings", "font", self.font)

        config.add_section("query_settings")
        config.set("query_settings", "query_data_radius", self.query_data_radius)
        config.set("query_settings", "query_phot_radius", self.query_phot_radius)
        config.set(
            "query_settings", "query_bulkphot_radius", self.query_bulkphot_radius
        )
        config.set(
            "query_settings", "query_lightcurve_radius", self.query_lightcurve_radius
        )
        config.set(
            "query_settings", "query_spectrum_radius", self.query_spectrum_radius
        )
        config.set("query_settings", "query_sed_radius", self.query_sed_radius)
        config.set(
            "query_settings", "query_reddening_radius", self.query_reddening_radius
        )
        config.set("query_settings", "query_image_size", self.query_image_size)
        config.set("query_settings", "query_image_overlays", self.query_image_overlays)
        config.set("query_settings", "query_image_band", self.query_image_band)
        config.set(
            "query_settings",
            "query_lightcurve_atlas_username",
            self.query_lightcurve_atlas_username,
        )
        config.set(
            "query_settings",
            "query_lightcurve_atlas_password",
            self.query_lightcurve_atlas_password,
        )

        config.add_section("image_overlay_settings")
        config.set("image_overlay_settings", "gaia_overlay_mag", self.gaia_overlay_mag)
        config.set(
            "image_overlay_settings", "galex_overlay_mag", self.galex_overlay_mag
        )
        config.set("image_overlay_settings", "wise_overlay_mag", self.wise_overlay_mag)
        config.set("image_overlay_settings", "sdss_overlay_mag", self.sdss_overlay_mag)
        config.set(
            "image_overlay_settings", "twomass_overlay_mag", self.twomass_overlay_mag
        )
        config.set(
            "image_overlay_settings",
            "skymapper_overlay_mag",
            self.skymapper_overlay_mag,
        )
        config.set(
            "image_overlay_settings",
            "panstarrs_overlay_mag",
            self.panstarrs_overlay_mag,
        )
        config.set(
            "image_overlay_settings",
            "overlay_piggyback_radius",
            self.overlay_piggyback_radius,
        )

        config.add_section("search_settings")
        config.set("search_settings", "search_radius", self.search_radius)

        config.add_section("datapage_settings")
        config.set(
            "datapage_settings",
            "datapage_search_button_radius",
            self.datapage_search_button_radius,
        )
        config.set(
            "datapage_settings",
            "datapage_datatable_radius",
            self.datapage_datatable_radius,
        )
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
        for key, val in vars(self).items():
            if key != "config_file":
                print(f"{key} = {val}")
