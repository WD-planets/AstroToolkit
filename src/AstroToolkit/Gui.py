"""
This module contains the ATK GUI, which allows for much of the package to be used through an interface.
"""

import os
import re
import sys

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (QApplication, QComboBox, QFileDialog, QLabel,
                             QLineEdit, QMainWindow, QPushButton)

from .Configuration.baseconfig import ConfigStruct
from .Data.dataquery import SurveyInfo
from .Tools import query, readdata

config = ConfigStruct()
config.read_config()

supported_surveys = {
    "data": SurveyInfo().list + ["other"],
    "phot": SurveyInfo().bulkphot_surveys,
    "reddening": SurveyInfo().reddening_surveys,
    "image": SurveyInfo().image_surveys,
    "lightcurve": SurveyInfo().lightcurve_surveys,
    "spectrum": SurveyInfo().spectrum_surveys,
}

grid_width = 150
grid_height = 50
margin = 5

newline = "\n"


class Window(QMainWindow):
    """
    :meta private:
    """

    def __init__(self):
        super().__init__()

        # Initialise initial values
        self.query_kind = "data"
        self.survey = "gaia"
        self.source = None
        self.sources = None
        self.pos = None
        self.username = None
        self.password = None
        self.data_extension = ".csv"
        self.data_ready = False

        self.setWindowTitle("ATK GUI")
        self.setGeometry(100, 100, 11 * grid_width, 20 * grid_height)
        self.MainUiComponents()
        self.QuerySettingsUiComponents()
        self.DataManipulationUiComponents()
        self.PlottingUiComponents()
        self.PlotOptionsUiComponents()
        self.set_additional_settings()
        self.set_plotoptions()
        self.update_button_positions()
        self.update_plotting_button_positions()

        self.show()

    # Sets position in grid space, and scales widgets to grid size
    def setPosition(self, widget, position):
        x, y = position

        widget.setGeometry(
            margin + x * grid_width, margin + y * grid_height, grid_width, grid_height
        )

    def getPosition(self, widget):
        return [
            int((widget.pos().x() - margin) / grid_width),
            int((widget.pos().y() - margin) / grid_height),
        ]

    def centerList(self, widget):
        widget.setEditable(True)
        line_edit = widget.lineEdit()
        line_edit.setAlignment(Qt.AlignCenter)
        line_edit.setReadOnly(True)

    def customLabel(self, text):
        label = QLabel(self)
        label.setText(text)
        header_font = QFont()
        header_font.setBold(True)
        label.setFont(header_font)
        label.setAlignment(Qt.AlignBottom | Qt.AlignCenter)
        label.setStyleSheet("padding :3px")

        return label

    def customLabel2(self, text):
        label = QLabel(self)
        label.setText(text)
        label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        label.setStyleSheet("padding :3px")

        return label

    def parseTextBox(self, textbox):
        arr = [a for a in re.split(r"(\s|\,)", textbox.text().strip()) if a]
        arr = [x for x in arr if x != " " and x != ","]
        return arr

    def MainUiComponents(self):
        """
        Main Query UI
        """

        # Labels
        self.query_kind_header = self.customLabel("Query Kind:")
        self.setPosition(self.query_kind_header, [0, 0])
        self.survey_header = self.customLabel("Survey:")
        self.setPosition(self.survey_header, [2, 0])

        # Query kind
        self.query_kind_list = QComboBox(self)
        self.query_kind_list.addItems(SurveyInfo().supported_query_kinds)
        self.setPosition(self.query_kind_list, [0, 1])
        self.query_kind_list.currentTextChanged.connect(self.on_change_linked)

        # Target Box
        self.target_box = QLineEdit(self, placeholderText="Target Source / Position")
        self.setPosition(self.target_box, [1, 1])

        # Survey Box
        self.survey_list = QComboBox(self)
        surveys = supported_surveys[self.query_kind]
        surveys = [survey for survey in surveys if survey != "gaia_lc"]
        self.survey_list.addItems(surveys)
        self.setPosition(self.survey_list, [2, 1])

        # Execute Query
        self.query_button = QPushButton(self, text="Execute Query")
        self.setPosition(self.query_button, [3, 1])

        # Vizier Catalogue
        self.vizier_catalogue = QLineEdit(self, placeholderText="Vizier Catalogue ID")
        self.setPosition(self.vizier_catalogue, [3, 1])
        self.vizier_catalogue.hide()

        # ATLAS username
        self.atlas_username = QLineEdit(self, placeholderText="ATLAS Username")
        self.setPosition(self.atlas_username, [3, 1])
        self.atlas_username.hide()

        # ATLAS password
        self.atlas_password = QLineEdit(self, placeholderText="ATLAS Password")
        self.setPosition(self.atlas_password, [4, 1])
        self.atlas_password.hide()

        # HRD sources
        self.hrd_sources = QLineEdit(self, placeholderText="Source(s)")
        self.setPosition(self.hrd_sources, [1, 1])
        self.hrd_sources.hide()

        # Read data
        self.readdata_button = QPushButton(self, text="Open File Browser")
        self.setPosition(self.readdata_button, [2, 4])

        self.readdata_header = self.customLabel(f"Read{newline}Local ATK File:")
        self.setPosition(self.readdata_header, [2, 3])

        # Signals
        self.query_kind_list.currentTextChanged.connect(self.on_change_linked)
        self.survey_list.currentTextChanged.connect(self.on_change_linked)

        self.query_kind_list.currentTextChanged.connect(self.on_query_kind_change)

        self.vizier_catalogue.textChanged.connect(self.on_vizier_catalogue_change)

        self.target_box.textChanged.connect(self.on_target_change)

        self.query_button.pressed.connect(self.execute_query)

        self.atlas_username.textChanged.connect(self.on_atlas_login_change)
        self.atlas_password.textChanged.connect(self.on_atlas_login_change)

        self.hrd_sources.textChanged.connect(self.on_hrd_sources_change)

        self.readdata_button.pressed.connect(self.get_file)

    # Read local files
    def get_file(self):
        path = os.getcwd()
        fname, _ = QFileDialog.getOpenFileName(
            self, "Open File", str(path), "ATK Data Files (*.csv *.fits)"
        )
        if fname:
            self.data = readdata(fname=fname)
            if hasattr(self.data, "subkind"):
                self.query_kind_list.setCurrentText(self.data.subkind)
            else:
                self.query_kind_list.setCurrentText(self.data.kind)

            if hasattr(self.data, "survey"):
                self.survey_list.setCurrentText(self.data.survey)
            if self.data.source:
                self.target_box.setText(str(self.data.source))
            else:
                self.target_box.setText(f"{self.pos[0]},{self.pos[1]}")
            self.on_data_found()

    # Triggers on both query kind and survey changes
    def on_change_linked(self):
        self.query_kind = self.query_kind_list.currentText()
        self.survey = self.survey_list.currentText()

        # Hide by default, shown if needed
        self.atlas_username.hide()
        self.atlas_password.hide()
        self.vizier_catalogue.hide()
        self.target_box.show()
        self.hrd_sources.hide()

        if self.query_kind == "data" and self.survey == "other":
            self.setPosition(self.query_button, [4, 1])
            self.vizier_catalogue.show()
        elif self.query_kind == "lightcurve" and self.survey == "atlas":
            self.setPosition(self.query_button, [5, 1])
            self.atlas_username.show()
            self.atlas_password.show()
        elif self.query_kind == "hrd":
            self.target_box.hide()
            self.hrd_sources.show()
        else:
            self.setPosition(self.query_button, [3, 1])

        if self.query_kind != "image":
            self.data_extension = ".csv"
        else:
            self.data_extension = ".fits"

        self.query_button.setStyleSheet("background-color : white")

    # handles hiding of survey button if not needed, and updates survey list
    def on_query_kind_change(self):
        if self.query_kind in ["bulkphot", "sed", "hrd"]:
            self.survey_list.hide()
            self.survey_header.hide()
            self.setPosition(self.query_button, [2, 1])
        else:
            surveys = supported_surveys[self.query_kind]
            if self.query_kind == "data":
                surveys = [survey for survey in surveys if survey != "gaia_lc"]
            self.survey_list.clear()
            self.survey_list.addItems(surveys)
            self.survey_list.show()
            self.survey_header.show()

    def on_vizier_catalogue_change(self):
        if not self.vizier_catalogue.isHidden():
            self.survey = self.vizier_catalogue.text()

    def on_target_change(self):
        target = self.target_box.text()

        # if target is a source
        if re.search("^[0-9]+$", target):
            self.source = int(target)
            self.pos = None
        else:
            target = self.parseTextBox(self.target_box)
            try:
                if len(target) == 2:
                    self.pos = [float(target[0]), float(target[1])]
                else:
                    self.pos = [0.0, 0.0]
            except:
                self.pos = [0.0, 0.0]
            self.source = None

    def on_hrd_sources_change(self):
        self.sources = self.parseTextBox(self.hrd_sources)

    def on_atlas_login_change(self):
        self.username = self.atlas_username.text()
        self.password = self.atlas_password.text()

    def execute_query(self):
        if self.query_kind in ["bulkphot", "sed", "hrd"]:
            self.survey = None
        if self.query_kind not in ["hrd"]:
            self.sources = None

        self.data = query(
            kind=self.query_kind,
            survey=self.survey,
            source=self.source,
            pos=self.pos,
            radius=self.radius,
            username=self.username,
            password=self.password,
            raw=self.return_raw,
            size=self.image_size,
            band=self.band,
            overlays=self.overlays,
            sources=self.sources,
        )

        self.on_data_found()

    def on_data_found(self):
        success = False
        if self.query_kind in ["lightcurve"]:
            for band in self.data.data:
                if band["mag"]:
                    success = True
        elif self.query_kind in ["bulkphot"]:
            for survey in self.data.data:
                if self.data.data[survey] is not None:
                    success = True
        else:
            if self.data.data:
                success = True

        if success:
            self.data_ready = True

            self.query_button.setStyleSheet("background-color : rgba(0,255,0,0.75)")
            self.savedata_fname_box.setPlaceholderText(self.data.dataname)
            if self.query_kind not in ["data", "reddening", "phot", "bulkphot"]:
                self.saveplot_fname_box.setPlaceholderText(self.data.plotname)
                self.set_plotoptions()
            for widget in self.data_manipulation_widgets:
                widget[0].show()
            if self.query_kind not in ["data", "phot", "bulkphot", "reddening"]:
                for widget in self.plotting_widgets:
                    widget[0].show()
        else:
            self.query_button.setStyleSheet("background-color: rgba(255,0,0,0.75)")

    # ADDITIONAL QUERY SETTINGS ----------------------------------------------------------------------------

    def QuerySettingsUiComponents(self):
        self.query_settings_header = self.customLabel(
            f"Additional{newline}Query Settings:"
        )
        self.setPosition(self.query_settings_header, [0, 3])

        self.additional_settings = {}
        self.additional_settings_headers = {}
        self.default_values = {}

        # radius
        query_radius = QLineEdit(self)
        query_radius_header = self.customLabel2(f"Query radius in{newline}arcseconds")
        self.additional_settings["query_radius"] = query_radius
        self.additional_settings_headers["query_radius"] = query_radius_header
        self.default_values["query_radius"] = getattr(
            config, f"query_{self.query_kind}_radius"
        )
        self.additional_settings["query_radius"].textChanged.connect(
            self.on_query_radius_change
        )
        self.radius = "config"

        # lightcurve raw
        return_raw = QComboBox(self)
        return_raw.addItems(["False", "True"])
        return_raw_header = self.customLabel2(
            f"Return raw data{newline}with no filtering"
        )
        self.additional_settings["return_raw"] = return_raw
        self.additional_settings_headers["return_raw"] = return_raw_header
        self.default_values["return_raw"] = False
        self.additional_settings["return_raw"].currentTextChanged.connect(
            self.on_return_raw_changed
        )
        self.return_raw = False

        # image size
        image_size = QLineEdit(self)
        image_size_header = self.customLabel2(f"Image size in{newline}arcseconds")
        self.additional_settings["image_size"] = image_size
        self.additional_settings_headers["image_size"] = image_size_header
        self.default_values["image_size"] = config.query_image_size
        self.additional_settings["image_size"].textChanged.connect(
            self.on_image_size_changed
        )
        self.image_size = "config"

        # image band
        image_band = QLineEdit(self)
        image_band_header = self.customLabel2(f"Image band as a{newline}single string")
        self.additional_settings["image_band"] = image_band
        self.additional_settings_headers["image_band"] = image_band_header
        self.default_values["image_band"] = config.query_image_band
        self.additional_settings["image_band"].textChanged.connect(
            self.on_image_band_changed
        )
        self.band = "config"

        # image overlays
        image_overlays = QLineEdit(self)
        image_overlays_header = self.customLabel2(
            f"Sets which overlay{newline}surveys to use."
        )
        self.additional_settings["image_overlays"] = image_overlays
        self.additional_settings_headers["image_overlays"] = image_overlays_header
        self.default_values["image_overlays"] = config.query_image_overlays
        self.additional_settings["image_overlays"].textChanged.connect(
            self.on_image_overlays_changed
        )
        self.overlays = "config"

        self.query_kind_list.currentTextChanged.connect(self.set_additional_settings)
        self.survey_list.currentTextChanged.connect(self.set_additional_settings)

    def on_query_radius_change(self):
        if self.additional_settings["query_radius"].text():
            self.radius = float(self.additional_settings["query_radius"].text())
        else:
            self.radius = "config"

    def on_return_raw_changed(self):
        if self.additional_settings["return_raw"].currentText().lower() == "true":
            self.return_raw = True
        else:
            self.return_Raw = False

    def on_image_size_changed(self):
        if self.additional_settings["image_size"].text():
            self.size = int(self.additional_settings["image_size"].text())
        else:
            self.size = "config"

    def on_image_band_changed(self):
        if self.additional_settings["image_band"].text():
            self.band = str(self.additional_settings["image_band"].text())
        else:
            self.band = "config"

    def on_image_overlays_changed(self):
        if self.additional_settings["image_overlays"].text():
            self.overlays = self.parseTextBox(
                self.additional_settings["image_overlays"]
            )
        else:
            self.overlays = "config"

    # triggers on query kind change and survey change
    def set_additional_settings(self):
        for key, val in self.additional_settings.items():
            if key not in ["return_raw"]:
                val.hide()
                val.clear()
            else:
                val.setCurrentText(str(self.default_values[key]))
                val.hide()
        for key, val in self.additional_settings_headers.items():
            val.hide()

        if self.query_kind == "data":
            params = ["query_radius"]
        elif self.query_kind == "phot":
            params = ["query_radius"]
        elif self.query_kind == "bulkphot":
            params = ["query_radius"]
        elif self.query_kind == "reddening":
            if self.survey == "stilism":
                params = []
            else:
                params = ["query_radius"]
        elif self.query_kind == "image":
            if self.survey in ["panstarrs", "skymapper"]:
                params = ["image_size", "image_band", "image_overlays"]
            else:
                params = ["image_size", "image_overlays"]
        elif self.query_kind == "lightcurve":
            params = ["query_radius", "return_raw"]
        elif self.query_kind == "hrd":
            params = []
        elif self.query_kind == "sed":
            params = ["query_radius"]
        elif self.query_kind == "spectrum":
            params = ["query_radius"]

        if len(params) == 0:
            self.query_settings_header.hide()
        else:
            self.query_settings_header.show()

        for index, key in enumerate(params):
            if key == "query_radius":
                self.default_values[key] = getattr(
                    config, f"query_{self.query_kind}_radius"
                )

            self.additional_settings[key].setPlaceholderText(
                str(self.default_values[key])
            )

            self.setPosition(self.additional_settings[key], [0, index + 4])
            self.setPosition(self.additional_settings_headers[key], [1, index + 4])
            self.additional_settings[key].show()
            self.additional_settings_headers[key].show()

    # DATA MANIPULATION ------------------------------------------------------------------------------------

    def DataManipulationUiComponents(self):
        self.data_manipulation_widgets = []

        self.showdata_button = QPushButton(self, text="Show Data")
        self.showdata_button.pressed.connect(self.on_showdata_button_pressed)
        self.data_manipulation_widgets.append([self.showdata_button, 2, 0])

        self.savedata_button = QPushButton(self, text="Save Data")
        self.savedata_button.pressed.connect(self.on_savedata_button_pressed)
        self.data_manipulation_widgets.append([self.savedata_button, 4, 0])

        self.showraw_box = QComboBox(self)
        self.showraw_box.addItems(["True", "False"])
        self.showraw_box.currentTextChanged.connect(self.on_showraw_box_change)
        self.showraw = False
        self.data_manipulation_widgets.append([self.showraw_box, 2, 3])

        self.showraw_header = self.customLabel(f"Additional{newline}Settings:")
        self.data_manipulation_widgets.append([self.showraw_header, 2, 2])

        self.showraw_header2 = self.customLabel2("Readable output")
        self.data_manipulation_widgets.append([self.showraw_header2, 3, 3])

        self.savedata_fname_box = QLineEdit(self)
        self.data_manipulation_widgets.append([self.savedata_fname_box, 4, 3])

        self.savedata_fname_header = self.customLabel(f"Additional{newline}Settings:")
        self.data_manipulation_widgets.append([self.savedata_fname_header, 4, 2])

        self.savedata_fname_header2 = self.customLabel2("File name")
        self.data_manipulation_widgets.append([self.savedata_fname_header2, 5, 3])

        self.query_kind_list.currentTextChanged.connect(self.update_button_positions)
        self.survey_list.currentTextChanged.connect(self.update_button_positions)
        self.savedata_fname_box.textChanged.connect(self.on_savedata_fname_changed)

    def update_button_positions(self):
        for widget in self.data_manipulation_widgets:
            widget[0].hide()

        # self.savedata_button.setText("Save Data")
        # self.savedata_button.setStyleSheet("background-color: white")

        query_button_position = self.getPosition(self.query_button)
        for widget in self.data_manipulation_widgets:
            self.setPosition(
                widget[0],
                [
                    query_button_position[0] + widget[1],
                    query_button_position[1] + widget[2],
                ],
            )

    def on_showdata_button_pressed(self):
        self.data.showdata(raw=self.showraw)

    def on_savedata_button_pressed(self):
        if self.savedata_fname_box.text():
            fname = self.savedata_fname_box.text()
            fname = self.data.savedata(name=fname)
        else:
            fname = self.data.savedata()

        if fname not in keep_files:
            keep_files.append(fname)

        self.savedata_button.setText("Data Saved")
        self.savedata_button.setStyleSheet("background-color: rgba(0,255,0,0.75)")

    def on_savedata_fname_changed(self):
        self.savedata_button.setText("Save Data")
        self.savedata_button.setStyleSheet("background-color: white")

        if self.savedata_fname_box.text():
            fname = self.savedata_fname_box.text()
            if not fname.endswith(self.data_extension):
                fname += self.data_extension
            if fname in keep_files:
                self.savedata_button.setText("Data Saved")
                self.savedata_button.setStyleSheet(
                    "background-color: rgba(0,255,0,0.75)"
                )
        else:
            fname = self.savedata_fname_box.placeholderText()
            if not fname.endswith(self.data_extension):
                fname += self.data_extension
            if fname in keep_files:
                self.savedata_button.setText("Data Saved")
                self.savedata_button.setStyleSheet(
                    "background-color: rgba(0,255,0,0.75)"
                )

    def on_showraw_box_change(self):
        if self.showraw_box.currentText().lower() == "true":
            self.showraw = False
        else:
            self.showraw = True

    # PLOTTING ---------------------------------------------------------------------------------------------

    def PlottingUiComponents(self):
        self.plotting_widgets = []

        self.showplot_button = QPushButton(self, text="Plot Data")
        self.showplot_button.pressed.connect(self.on_showplot_button_pressed)
        self.plotting_widgets.append([self.showplot_button, 0, 7])

        self.saveplot_button = QPushButton(self, text="Save Plot")
        self.saveplot_button.pressed.connect(self.on_saveplot_button_pressed)
        self.plotting_widgets.append([self.saveplot_button, 2, 7])

        self.saveplot_fname_box = QLineEdit(self)
        self.plotting_widgets.append([self.saveplot_fname_box, 2, 10])
        self.saveplot_fname_box.textChanged.connect(self.on_saveplot_fname_changed)

        self.saveplot_fname_header = self.customLabel(f"Additional{newline}Settings:")
        self.plotting_widgets.append([self.saveplot_fname_header, 2, 9])

        self.saveplot_fname_header2 = self.customLabel2("File name")
        self.plotting_widgets.append([self.saveplot_fname_header2, 3, 10])

        self.query_kind_list.currentTextChanged.connect(
            self.update_plotting_button_positions
        )
        self.survey_list.currentTextChanged.connect(
            self.update_plotting_button_positions
        )

    def do_plotting(self):
        if self.query_kind == "lightcurve":
            plot_kind = self.plotoptions["kind"].currentText()
            timeformat = self.plotoptions["lightcurve_timeformat"].currentText()
            if plot_kind == "lightcurve":
                colours = self.parseTextBox(self.plotoptions["lightcurve_colours"])
                bands = self.parseTextBox(self.plotoptions["lightcurve_bands"])
                if self.plotoptions["lightcurve_sigmaclip"].text():
                    sigmaclip = float(self.plotoptions["lightcurve_sigmaclip"].text())
                    self.data.sigmaclip(sigmaclip)
                if self.plotoptions["lightcurve_bin"].text():
                    bin_text = self.plotoptions["lightcurve_bin"].text()
                    print(bin_text)
                    if (
                        bin_text.endswith("d")
                        or bin_text.endswith("h")
                        or bin_text.endswith("m")
                    ):
                        self.data.bin(binsize=bin_text)
                    else:
                        self.data.bin(bins=int(bin_text))
                self.data.plot(
                    kind=plot_kind,
                    colours=colours,
                    timeformat=timeformat,
                    bands=bands,
                )
            elif plot_kind == "powspec":
                method = self.plotoptions["powspec_method"].currentText()
                self.data.plot(kind="powspec", method=method)
            elif plot_kind == "phasefold":
                if self.plotoptions["phasefold_freq"].text():
                    freq = float(self.plotoptions["phasefold_freq"].text())
                else:
                    freq = None
                if self.plotoptions["phasefold_bins"].text():
                    bins = int(self.plotoptions["phasefold_bins"].text())
                else:
                    bins = None
                if (
                    self.plotoptions["phasefold_foverlay"].currentText().lower()
                    == "true"
                ):
                    foverlay = True
                else:
                    foverlay = False
                if self.plotoptions["phasefold_repeat"].text():
                    repeat = int(self.plotoptions["phasefold_repeat"].text())
                else:
                    repeat = self.plotoptions_defaults["phasefold_repeat"]
                if self.plotoptions["phasefold_shift"].text():
                    shift = float(self.plotoptions["phasefold_shift"].text())
                else:
                    shift = self.plotoptions_defaults["phasefold_shift"]
                self.data.plot(
                    kind="phasefold",
                    freq=freq,
                    bins=bins,
                    foverlay=foverlay,
                    repeat=repeat,
                    shift=shift,
                )
        elif self.query_kind == "sed":
            if self.plotoptions["sed_spec"].currentText() != "---":
                self.data.plot(
                    spectrum_overlay=True,
                    survey=self.plotoptions["sed_spec"].currentText(),
                )
            else:
                self.data.plot()
        else:
            self.data.plot()

    def on_showplot_button_pressed(self):
        self.do_plotting()
        fname = self.data.showplot()

        if fname not in files:
            files.append(fname)

    def on_saveplot_button_pressed(self):
        self.do_plotting()
        fname = self.data.saveplot()

        if fname not in keep_files:
            keep_files.append(fname)

        self.saveplot_button.setText("Plot Saved")
        self.saveplot_button.setStyleSheet("background-color: rgba(0,255,0,0.75)")

    def on_saveplot_fname_changed(self):
        self.saveplot_button.setText("Save Plot")
        self.saveplot_button.setStyleSheet("background-color: white")

        if self.saveplot_fname_box.text():
            fname = self.saveplot_fname_box.text()
            if not fname.endswith(".html"):
                fname += ".html"
            if fname in keep_files:
                self.saveplot_button.setText("Data Saved")
                self.saveplot_button.setStyleSheet(
                    "background-color: rgba(0,255,0,0.75)"
                )
        else:
            fname = self.saveplot_fname_box.placeholderText()
            if not fname.endswith(".html"):
                fname += ".html"
            if fname in keep_files:
                self.saveplot_button.setText("Data Saved")
                self.saveplot_button.setStyleSheet(
                    "background-color: rgba(0,255,0,0.75)"
                )

    def update_plotting_button_positions(self):
        for widget in self.plotting_widgets:
            widget[0].hide()

        query_kind_position = self.getPosition(self.query_kind_list)
        for widget in self.plotting_widgets:
            self.setPosition(
                widget[0],
                [
                    query_kind_position[0] + widget[1],
                    query_kind_position[1] + widget[2],
                ],
            )

    # PLOTTING OPTIONS -------------------------------------------------------------------------------------

    def PlotOptionsUiComponents(self):
        self.plotoptions_header = self.customLabel(
            f"Additional{newline}Plotting Options:"
        )
        self.setPosition(self.plotoptions_header, [0, 10])

        self.plotoptions = {}
        self.plotoptions_headers = {}
        self.plotoptions_defaults = {}

        # Lightcurve - kind
        plot_kind = QComboBox(self)
        plot_kind.addItems(["lightcurve", "powspec", "phasefold"])
        self.plotoptions["kind"] = plot_kind
        self.plotoptions_headers["kind"] = self.customLabel2("Plot type")
        self.plotoptions_defaults["kind"] = "lightcurve"

        # Lightcurve - bands
        self.plotoptions["lightcurve_bands"] = QLineEdit(self)
        self.plotoptions_headers["lightcurve_bands"] = self.customLabel2(
            "Light curve bands"
        )

        # Lightcurve - colours
        self.plotoptions["lightcurve_colours"] = QLineEdit(self)
        self.plotoptions_headers["lightcurve_colours"] = self.customLabel2(
            "Light curve colours"
        )

        # Lightcurve - timeformat
        self.plotoptions["lightcurve_timeformat"] = QComboBox(self)
        self.plotoptions["lightcurve_timeformat"].addItems(["reduced", "original"])
        self.plotoptions_headers["lightcurve_timeformat"] = self.customLabel2(
            "Time format"
        )
        self.plotoptions_defaults["lightcurve_timeformat"] = "reduced"

        # Lightcurve - sigmaclip()
        self.plotoptions["lightcurve_sigmaclip"] = QLineEdit(self)
        self.plotoptions_headers["lightcurve_sigmaclip"] = self.customLabel2(
            "Sigma clip"
        )
        self.plotoptions_defaults["lightcurve_sigmaclip"] = "---"

        # Lightcurve - bin()
        self.plotoptions["lightcurve_bin"] = QLineEdit(self)
        self.plotoptions_headers["lightcurve_bin"] = self.customLabel2(
            f"Bins or binsize{newline}e.g. 30 or 1h"
        )
        self.plotoptions_defaults["lightcurve_bin"] = "---"

        # Powspec - method
        self.plotoptions["powspec_method"] = QComboBox(self)
        self.plotoptions["powspec_method"].addItems(
            ["ls", "amhw", "pspw", "atrw", "aovw", "f_mw", "lomw"]
        )
        self.plotoptions_headers["powspec_method"] = self.customLabel2(
            "Analysis Method"
        )
        self.plotoptions_defaults["powspec_method"] = "ls"

        # Phasefold - freq
        self.plotoptions["phasefold_freq"] = QLineEdit(self)
        self.plotoptions_headers["phasefold_freq"] = self.customLabel2("Fold frequency")
        self.plotoptions_defaults["phasefold_freq"] = "Calculate"

        # Phasefold - bins
        self.plotoptions["phasefold_bins"] = QLineEdit(self)
        self.plotoptions_headers["phasefold_bins"] = self.customLabel2("Number of bins")
        self.plotoptions_defaults["phasefold_bins"] = "---"

        # Phasefold - foverlay
        self.plotoptions["phasefold_foverlay"] = QComboBox(self)
        self.plotoptions["phasefold_foverlay"].addItems(["True", "False"])
        self.plotoptions_headers["phasefold_foverlay"] = self.customLabel2(
            "Sine wave overlay"
        )
        self.plotoptions_defaults["phasefold_foverlay"] = True

        # Phasefold - repeat
        self.plotoptions["phasefold_repeat"] = QLineEdit(self)
        self.plotoptions_headers["phasefold_repeat"] = self.customLabel2(
            "Number of repetitions"
        )
        self.plotoptions_defaults["phasefold_repeat"] = 2

        # Phasefold - shift587316166180416640
        self.plotoptions["phasefold_shift"] = QLineEdit(self)
        self.plotoptions_headers["phasefold_shift"] = self.customLabel2(
            "Shift in phase"
        )
        self.plotoptions_defaults["phasefold_shift"] = 0

        # SED - spectrum_overlay
        self.plotoptions["sed_spec"] = QComboBox(self)
        self.plotoptions["sed_spec"].addItems(["---"] + SurveyInfo().spectrum_surveys)
        self.plotoptions_headers["sed_spec"] = self.customLabel2(
            "Spectrum overlay survey"
        )
        self.plotoptions_defaults["sed_spec"] = "---"

        self.survey_list.currentTextChanged.connect(self.disable_widgets)
        self.query_kind_list.currentTextChanged.connect(self.set_plotoptions)
        self.survey_list.currentTextChanged.connect(self.set_plotoptions)
        self.plotoptions["kind"].currentTextChanged.connect(self.set_plotoptions, True)

    def set_plotoptions(self, skip_reset=False):
        self.saveplot_button.setStyleSheet("background-color: white")
        self.saveplot_button.setText("Save Plot")

        for key, val in self.plotoptions.items():
            # put names of QComboBoxes here so that they don't get cleared
            if key not in [
                "kind",
                "lightcurve_timeformat",
                "powspec_method",
                "phasefold_foverlay",
                "sed_spec",
            ]:
                val.hide()
                val.clear()
            else:
                # this also gets called on change of plot kind, so need to not reset this to default in this case
                if not skip_reset:
                    val.setCurrentText(str(self.plotoptions_defaults[key]))
                val.hide()
        for key, val in self.plotoptions_headers.items():
            val.hide()

        if self.data_ready:
            self.plotoptions_header.show()

            if self.query_kind == "lightcurve":
                # place defaults here that need to update per-query
                lightcurve_bands = ""
                for band in SurveyInfo().lightcurve_bands[self.survey]:
                    lightcurve_bands += band + ", "
                lightcurve_bands = lightcurve_bands[:-2]
                self.plotoptions_defaults["lightcurve_bands"] = lightcurve_bands

                lightcurve_colours = (
                    "black, " * len(SurveyInfo().lightcurve_bands[self.survey])
                )[:-2]
                self.plotoptions_defaults["lightcurve_colours"] = lightcurve_colours

                if self.plotoptions["kind"].currentText() == "lightcurve":
                    params = [
                        "kind",
                        "lightcurve_bands",
                        "lightcurve_colours",
                        "lightcurve_timeformat",
                        "lightcurve_sigmaclip",
                        "lightcurve_bin",
                    ]
                elif self.plotoptions["kind"].currentText() == "powspec":
                    params = ["kind", "powspec_method"]
                elif self.plotoptions["kind"].currentText() == "phasefold":
                    params = [
                        "kind",
                        "phasefold_freq",
                        "phasefold_bins",
                        "phasefold_foverlay",
                        "phasefold_repeat",
                        "phasefold_shift",
                    ]
            elif self.query_kind == "image":
                params = []
            elif self.query_kind == "sed":
                params = ["sed_spec"]
            elif self.query_kind == "spectrum":
                params = []
            elif self.query_kind == "hrd":
                params = []

            # Don't have plotting, but included here so that it doesn't break if these kinds are selected
            elif self.query_kind == "data":
                params = []
            elif self.query_kind == "phot":
                params = []
            elif self.query_kind == "bulkphot":
                params = []
            elif self.query_kind == "reddening":
                params = []

            if len(params) == 0:
                self.plotoptions_header.hide()
            else:
                self.plotoptions_header.show()

            for index, key in enumerate(params):
                self.plotoptions[key].setPlaceholderText(
                    str(self.plotoptions_defaults[key])
                )

                self.setPosition(self.plotoptions[key], [0, index + 11])
                self.setPosition(self.plotoptions_headers[key], [1, index + 11])
                self.plotoptions[key].show()
                self.plotoptions_headers[key].show()
        else:
            self.plotoptions_header.hide()

    def disable_widgets(self):
        self.data_ready = False


files = []
keep_files = []


def appExec():
    """
    :meta private:
    """
    app = QApplication([])
    window = Window()
    app.exec_()

    import os

    for file in files:
        if file not in keep_files:
            os.remove(file)


def openGUI():
    """
    Opens the ATK GUI. This function can also be called directly from the :ref:`Command-Line`.

    :return: None

    |

    """
    sys.exit(appExec())
