def check_targeting(input):
    if "pos" in input and "source" in input:
        if input["pos"] and input["source"]:
            raise ValueError(
                "Simultaneous source and pos input detected. Only one may be used."
            )
        elif not input["pos"] and not input["source"]:
            raise ValueError("Source or pos input required.")
    return input


def check_pos(pos):
    if not isinstance(pos, list):
        raise ValueError("Invalid pos input.")
    for coord in pos:
        if not isinstance(coord, int) and not isinstance(coord, float):
            try:
                coord = float(coord)
            except:
                raise ValueError("Invalid pos input.")
    return pos


def check_source(source):
    if not isinstance(source, int):
        try:
            source = int(source)
        except:
            raise ValueError("Invalid source input.")
    return source


def check_radius(radius):
    if not isinstance(radius, int) and not isinstance(radius, float):
        try:
            radius = float(radius)
        except:
            raise ValueError("Invalid radius input.")
    return radius


def check_size(size, survey):
    if not isinstance(size, int):
        try:
            size = int(size)
        except:
            raise ValueError("Invalid size input")
    if survey == "panstarrs" and size > 1500:
        raise ValueError('Size too large. Maximum supported by panstarrs is 1500".')
    elif survey == "skymapper" and size > 600:
        raise ValueError('Size too large. Maximum supported by skymapper is 600".')
    elif survey == "dss" and size > 7200:
        raise ValueError('Size too large. Maximum supported by dss is 7200".')
    return size


def check_band(band, survey):
    import re

    if survey == "panstarrs" and not re.match("^[grizy]+$", band):
        raise ValueError("Invalid panstarrs bands. Supported bands are [g,r,i,z,y].")
    elif survey == "skymapper":
        if not re.match("^[grizuv]+$", band):
            raise ValueError(
                "Invalid skymapper bands. Supported bands are [g,r,i,z,u,v]."
            )
        else:
            band = list(band)
            temp_string = ""
            for i in range(0, len(band)):
                temp_string += band[i] + ","
            band = temp_string[:-1]
    elif survey == "dss" and band != "g":
        raise ValueError("Invalid dss band. Only band supported by dss is g.")
    return band


def check_overlays(overlays):
    from ..Data.dataquery import SurveyInfo

    supported_overlays = SurveyInfo().supported_overlays

    if isinstance(overlays, list):
        for i, val in enumerate(overlays):
            overlays[i] = str(val).lower()
            if overlays[i] not in supported_overlays:
                raise ValueError(f"Unsupported overlay '{val}'.")
    elif isinstance(overlays, dict):
        for key, val in overlays.items():
            if key not in supported_overlays:
                raise ValueError(f"Unsupported overlay: '{val}'.")

    return overlays


def check_username(username, survey):
    if survey != "atlas":
        raise ValueError("Username only required in atlas queries.")
    if not isinstance(username, str):
        try:
            username = str(username)
        except:
            raise ValueError("Invalid atlas username.")
    return username


def check_password(password, survey):
    if survey != "atlas":
        raise ValueError("Password only required in atlas queries.")
    if not isinstance(password, str):
        try:
            password = str(password)
        except:
            raise ValueError("Invalid atlas password")
    return password


def check_time(time):
    if not isinstance(time, list):
        for element in time:
            if not isinstance(element, int):
                raise ValueError("Invalid time format.")
    return time


def check_value(value):
    if not isinstance(value, int) and not isinstance(value, float):
        try:
            value = float(value)
        except:
            raise ValueError(f"Invalid input. Expected int/float, got {type(value)}.")
    return value


def check_string(string):
    if not isinstance(string, str):
        try:
            string = str(string)
        except:
            raise ValueError(f"Invalid input. Expected str, got {type(string)}.")
    return string


def check_selection(selection):
    from ..Data.dataquery import SurveyInfo

    supported_surveys = SurveyInfo().list

    if not isinstance(selection, dict):
        raise ValueError(
            f"Invalid selection type. Expected dict, got {type(selection)}."
        )

    for survey, data in selection.items():
        if data != "default":
            if "parameters" not in data or "errors" not in data or "notes" not in data:
                raise ValueError("Input datatable selection missing required keys.")
            if survey not in supported_surveys:
                if "values" not in data:
                    raise ValueError("Input datatable selection missing required keys.")
    return selection


def check_dimensions(dimensions):
    if not isinstance(dimensions, dict):
        raise ValueError(f"Invalid input. Expected dict, got {type(dimensions)}.")

    if "width" not in dimensions or "height" not in dimensions:
        raise ValueError("Required keys for dimensions input: 'height', 'width'.")

    for key, val in dimensions.items():
        if not isinstance(val, int):
            try:
                dimensions[key] = int(val)
            except:
                raise ValueError("Dimensions must be integers.")
    return dimensions


def check_plots(plots):
    from bokeh.models.layouts import Column as bokeh_column
    from bokeh.models.layouts import Row as bokeh_row
    from bokeh.models.widgets.tables import DataTable as bokeh_datatable
    from bokeh.plotting._figure import figure as bokeh_figure

    from AstroToolkit.Data.hrdquery import HrdStruct
    from AstroToolkit.Data.imagequery import ImageStruct
    from AstroToolkit.Data.lightcurvequery import LightcurveStruct
    from AstroToolkit.Data.sedquery import SedStruct
    from AstroToolkit.Data.spectrumquery import SpectrumStruct

    for plot_info in plots:
        if (
            "name" not in plot_info
            or "width" not in plot_info
            or "height" not in plot_info
            or "figure" not in plot_info
        ):
            raise ValueError(
                "Required keys for grid entries: 'name', 'width', 'height', 'figure'."
            )
        if not isinstance(plot_info["name"], str):
            try:
                plot_info["name"] = str(plot_info["name"])
            except:
                raise ValueError("plot name must be string.")
        if not isinstance(plot_info["width"], int):
            try:
                plot_info["width"] = int(plot_info["width"])
            except:
                raise ValueError("width must int.")
        if not isinstance(plot_info["height"], int):
            try:
                plot_info["height"] = int(plot_info["height"])
            except:
                raise ValueError("height must be int.")
        if plot_info["figure"] and not isinstance(
            plot_info["figure"],
            (
                bokeh_figure,
                bokeh_row,
                bokeh_column,
                bokeh_datatable,
                ImageStruct,
                HrdStruct,
                LightcurveStruct,
                SedStruct,
                SpectrumStruct,
            ),
        ):
            raise ValueError(
                f"Note: Unsupported figure passed to gridsetup. Supported types are: ATK structure, bokeh row/column/figure/datatable or None, got {type(plot_info['figure'])}."
            )
    return plots


class InputObject(object):
    def __init__(self, input):
        self.input = input


class ValidateDataInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    self.input[key] = check_string(val)
                elif key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidateSpectrumInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    from ..Data.dataquery import SurveyInfo

                    supported_surveys = SurveyInfo().spectrum_surveys

                    if val not in supported_surveys:
                        raise ValueError("Invalid survey input.")
                elif key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidateImageInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    from ..Data.dataquery import SurveyInfo

                    supported_surveys = SurveyInfo().image_surveys

                    if val not in supported_surveys:
                        raise ValueError("Invalid survey input.")
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "size":
                    self.input[key] = check_size(val, self.input["survey"])
                elif key == "band":
                    self.input[key] = check_band(val, self.input["survey"])
                elif key == "overlays":
                    self.input[key] = check_overlays(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidateLightcurveInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    from ..Data.dataquery import SurveyInfo

                    supported_surveys = SurveyInfo().lightcurve_surveys

                    if val not in supported_surveys:
                        raise ValueError("Invalid survey input.")
                elif key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "username":
                    self.input[key] = check_username(val, self.input["survey"])
                elif key == "password":
                    self.input[key] = check_password(val, self.input["survey"])

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidatePhotInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    from ..Data.dataquery import SurveyInfo

                    supported_surveys = SurveyInfo().bulkphot_surveys

                    if val not in supported_surveys:
                        raise ValueError("Invalid survey input.")
                elif key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidateBulkphotInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateSedInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateReddeningInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "survey":
                    from ..Data.dataquery import SurveyInfo

                    supported_surveys = SurveyInfo().reddening_surveys

                    if val not in supported_surveys:
                        raise ValueError("Invalid survey input.")
                elif key == "radius":
                    if self.input["survey"] == "stilism":
                        raise ValueError(
                            "Radius not needed for stilism reddening query."
                        )
                    self.input[key] = check_radius(val)
                elif key == "pos":
                    if self.input["survey"] == "stilism":
                        raise ValueError("Stilism query does not support pos input.")
                    self.input[key] = check_pos(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        if not self.input["survey"]:
            raise ValueError("Survey required for query.")

        return self.input


class ValidateHrdInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "sources":
                    if isinstance(val, list):
                        for i, source in enumerate(val):
                            val[i] = check_source(source)
                        self.input[key] = val
                    else:
                        self.input[key] = check_source(val)

        if not self.input["sources"]:
            raise ValueError("Sources required for hrd query.")

        return self.input


class ValidateCorrectpmInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            from ..Data.dataquery import SurveyInfo

            supported_surveys = list(SurveyInfo().times.keys())

            if val and val != "config":
                if key == "input_survey":
                    if val not in supported_surveys:
                        raise ValueError("Invalid input_survey.")
                elif key == "target_survey":
                    if val not in supported_surveys:
                        raise ValueError("Invalid target_survey")
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "input_time":
                    self.input[key] = check_time(val)
                elif key == "target_time":
                    self.input[key] = check_time(val)
                elif key == "pmra":
                    self.input[key] = check_value(val)
                elif key == "pmdec":
                    self.input[key] = check_value(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateReaddataInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "fname":
                    self.input[key] = check_string(val)

        return self.input


class ValidateSearchInput(InputObject):
    def check_input(self):
        check_targeting(self.input)

        for key, val in self.input.items():
            if val and val != "config":
                if key == "kind":
                    if val not in ["simbad", "vizier"]:
                        raise ValueError("Invalid search kind.")
                elif key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateReadfitsInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "fname":
                    self.input[key] = check_string(val)
                elif key == "columns":
                    if isinstance(val, list):
                        for i, col in enumerate(val):
                            val[i] = check_string(col)
                        self.input[key] = val
                    else:
                        self.input[key] = check_string(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateDeg2hmsInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "pos":
                    self.input[key] = check_pos(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateHms2degInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "pos":
                    if not isinstance(val, str):
                        raise ValueError("pos must be str in hms2deg.")
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateButtonsInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "grid_size":
                    if not isinstance(val, int):
                        try:
                            self.input[key] = int(val)
                        except:
                            raise ValueError("grid_size must be int.")
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateDatatableInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "radius":
                    self.input[key] = check_radius(val)
                elif key == "source":
                    self.input[key] = check_source(val)
                elif key == "pos":
                    self.input[key] = check_pos(val)
                elif key == "selection":
                    check_selection(val)
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateEditInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "key":
                    from ..Configuration.baseconfig import ConfigStruct

                    config = ConfigStruct()
                    config.set_default_config()
                    accepted_keys = config.supported_keys

                    if val not in accepted_keys:
                        raise ValueError("Invalid config key.")
                elif key == "value":
                    pass
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        return self.input


class ValidateGridsetupInput(InputObject):
    def check_input(self):
        for key, val in self.input.items():
            if val and val != "config":
                if key == "dimensions":
                    self.input[key] = check_dimensions(val)
                elif key == "plots":
                    self.input[key] = check_plots(val)
                elif key == "grid_size":
                    if not isinstance(val, int):
                        try:
                            self.input[key] = int(val)
                        except:
                            raise ValueError("grid_size must be int.")
                else:
                    raise ValueError(f"Invalid input parameter '{key}'.")

        unit_area = 0
        for plot_info in self.input["plots"]:
            unit_area += plot_info["width"] * plot_info["height"]
        if (
            unit_area
            < self.input["dimensions"]["width"] * self.input["dimensions"]["height"]
        ):
            raise ValueError(
                "Given dimensions must be filled with figures. Pass entries with 'figure':None to fill empty space."
            )
        elif (
            unit_area
            > self.input["dimensions"]["width"] * self.input["dimensions"]["height"]
        ):
            raise ValueError(
                "Total area of elements is larger than the given dimensions."
            )

        return self.input


def check_inputs(input, kind):
    input_object = globals()[f"Validate{kind.capitalize()}Input"](input)
    return input_object.check_input()
