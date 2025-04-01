import configparser
import os

from ..PackageInfo import SurveyInfo


class EpochStruct(object):
    def __init__(self):
        from importlib_resources import files

        self.epoch_file = files("AstroToolkit.Configuration").joinpath("ATKEpochs.ini")
        if not os.path.isfile(self.epoch_file):
            print("No ATKEpochs.ini found. Generating one with default values...")
            self.default_setup()

    def default_setup(self):
        Epochs = configparser.ConfigParser()
        surveyInfo = SurveyInfo()

        default_data_surveys = surveyInfo.defaultDataSurveys
        default_data_surveys.remove("gaia_lc")
        lightcurve_surveys = surveyInfo.defaultLightcurveSurveys
        spectrum_surveys = surveyInfo.defaultSpectrumSurveys

        Epochs.add_section("default_data_surveys")
        Epochs.add_section("additional_data_surveys")
        Epochs.add_section("lightcurve_surveys")
        Epochs.add_section("spectrum_surveys")

        survey_epochs = surveyInfo.defaultSurveyTimes
        for survey in default_data_surveys:
            Epochs.set("default_data_surveys", survey, f"{survey_epochs[survey][0]},{survey_epochs[survey][1]}")
        for survey in lightcurve_surveys:
            Epochs.set("lightcurve_surveys", survey, f"{survey_epochs[survey][0]},{survey_epochs[survey][1]}")
        for survey in spectrum_surveys:
            Epochs.set("spectrum_surveys", survey, f"{survey_epochs[survey][0]},{survey_epochs[survey][1]}")

        with open(self.epoch_file, "w") as file:
            Epochs.write(file)

    def write_epochs(self):
        Epochs = configparser.ConfigParser()

        for section, entries in self.epochs.items():
            Epochs.add_section(section)
            for survey, epoch in entries.items():
                Epochs.set(section, survey, f"{epoch[0]},{epoch[1]}")

        with open(self.epoch_file, "w") as file:
            Epochs.write(file)

    def get_epochs(self):
        Epochs = configparser.ConfigParser()
        Epochs.read(self.epoch_file)

        sections = []
        sections.append(Epochs["default_data_surveys"])
        sections.append(Epochs["additional_data_surveys"])
        sections.append(Epochs["lightcurve_surveys"])
        sections.append(Epochs["spectrum_surveys"])

        section_strings = ["default_data_surveys", "additional_data_surveys", "lightcurve_surveys", "spectrum_surveys"]

        self.structured_out = {}
        self.epochs = {}
        for section, section_str in zip(sections, section_strings):
            self.structured_out[section_str] = []
            self.epochs[section_str] = {}
            for key, val in section.items():
                self.structured_out[section_str].append({key: val})
                epoch_list = val.split(",")
                self.epochs[section_str][key] = epoch_list

        return self.epochs

    @property
    def epoch_list(self):
        epoch_dict = self.get_epochs()

        epochs = {}
        for label, section in epoch_dict.items():
            for survey, epoch in section.items():
                if label == "lightcurve_surveys" and survey == "gaia":
                    epochs["gaia_lc"] = epoch
                else:
                    epochs[survey] = epoch

        for survey in epochs:
            epochs[survey] = [int(x) for x in epochs[survey]]

        return epochs

    def output_epochs(self):
        self.get_epochs()
        for label, section in self.structured_out.items():
            print(f"[{label}]")
            for entry in section:
                for key, val in entry.items():
                    print(f"{key} = {val}")
            print()

    def delete_epoch(self, survey):
        self.get_epochs()

        if survey in self.epochs["additional_data_surveys"]:
            del self.epochs["additional_data_surveys"][survey]
        else:
            raise ValueError(f"Could not find existing epoch definition '{survey}' in [additional_data_surveys].")

        self.write_epochs()

    def set_epoch(self, section, survey, epoch):
        self.get_epochs()
        if section not in ["default_data_surveys", "additional_data_surveys", "lightcurve_surveys", "spectrum_surveys"]:
            raise ValueError("Invalid section.")

        if survey in self.epochs[section]:
            print(f"Overwriting existing epoch for default survey '{survey}' in [{section}].")
        else:
            if section != "additional_data_surveys":
                raise ValueError(f"Cannot add new epoch definition in [{section}].")
            else:
                print(f"Adding new epoch definition for alias '{survey}' in [{section}].")

        self.epochs[section][survey] = epoch

        self.write_epochs()
