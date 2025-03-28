import configparser
import os

from ..PackageInfo import SurveyInfo

surveyInfo = SurveyInfo()


class CatalogueStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.catalogue_file = files("AstroToolkit.Configuration").joinpath(
            "ATKAliases.ini"
        )
        if not os.path.isfile(self.catalogue_file):
            print("No ATKAliases.ini found. Generating one with default values...")
            self.default_setup()

    def default_setup(self):
        Catalogues = configparser.ConfigParser()

        default_catalogues = surveyInfo.getCatalogueDict
        del default_catalogues["gaia_lc"]
        Catalogues.add_section("default_catalogues")
        Catalogues.add_section("additional_catalogues")
        for key, val in default_catalogues.items():
            Catalogues.set("default_catalogues", key, val)

        with open(self.catalogue_file, "w") as file:
            Catalogues.write(file)

    def get_catalogues(self):
        Catalogues = configparser.ConfigParser()
        Catalogues.read(self.catalogue_file)

        sections = []
        sections.append(Catalogues["default_catalogues"])
        sections.append(Catalogues["additional_catalogues"])

        sections_strings = ["default_catalogues", "additional_catalogues"]

        self.structured_out = {}
        for section, section_str in zip(sections, sections_strings):
            self.structured_out[section_str] = []
            for key, val in section.items():
                self.structured_out[section_str].append({key: val})
                setattr(self, key, val)

    def get_catalogue_list(self):
        self.get_catalogues()
        surveys = {}
        for survey, id in vars(self).items():
            if survey not in ["catalogue_file", "structured_out"]:
                surveys[survey] = id
        return surveys

    def get_alias_list(self):
        return [
            x for x in self.get_catalogue_list() if x not in surveyInfo.dataSurveyInfo
        ]

    def write_catalogues(self):
        Catalogues = configparser.ConfigParser()

        labels, ids = [], []
        for key, val in vars(self).items():
            if key not in ["catalogue_file", "structured_out"]:
                val = str(val)
                labels.append(key)
                ids.append(val)

        catalogues = surveyInfo.getCatalogueDict

        Catalogues.add_section("default_catalogues")
        Catalogues.add_section("additional_catalogues")
        for key, val in zip(labels, ids):
            if key in catalogues:
                Catalogues.set("default_catalogues", key, val)
            else:
                Catalogues.set("additional_catalogues", key, val)

        with open(self.catalogue_file, "w") as file:
            Catalogues.write(file)

    def delete_catalogue(self, key):
        self.get_catalogues()

        if hasattr(self, key):
            delattr(self, key)
        else:
            raise ValueError(f"Could not find Alias with label {key}.")

        self.write_catalogues()

    def add_catalogue(self, key, value):
        default_catalogues = surveyInfo.getCatalogueDict

        self.get_catalogues()
        key, value = str(key), str(value)
        if key in default_catalogues:
            raise ValueError(f"Cannot override default alias '{key}'.")
        if hasattr(self, key):
            print(
                f"Note: {key} alias was already defined, and has hence been overwritten."
            )
        setattr(self, key, value)
        self.write_catalogues()

    def output_catalogues(self):
        self.get_catalogues()
        for label, section in self.structured_out.items():
            print(f"[{label}]")
            for entry in section:
                for key, val in entry.items():
                    print(f"{key} = {val}")
            print()
