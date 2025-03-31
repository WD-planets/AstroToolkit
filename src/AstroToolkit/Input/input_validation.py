import re

from ..Configuration.overlays import OverlayStruct
from ..PackageInfo import SurveyInfo, ToolInfo

surveyInfo = SurveyInfo()
toolInfo = ToolInfo()

DO_NOT_LOWER = ["fname"]


def check_type(name, value, target_type):
    basic_checks_passed = False

    if not target_type:
        return value

    # check bools separately, as these should not be corrected due to ambiguity
    if target_type is bool and type(value) is not bool:
        raise ValueError(f"Invalid [{name}] input. Expected {target_type}, got {type(value)}.")

    # general type check, this already validates e.g. lists - further type checking happens below
    if not isinstance(value, target_type):
        if target_type is list:
            try:
                value = [value]
                basic_checks_passed = True
            except:
                basic_checks_passed = False
        else:
            try:
                value = target_type(value)
                basic_checks_passed = True
            except:
                basic_checks_passed = False
    else:
        basic_checks_passed = True

    if not basic_checks_passed:
        raise ValueError(f"Invalid [{name}] input. Expected {target_type}, got {type(value)}.")

    if target_type is str and name not in DO_NOT_LOWER:
        value = value.lower()

    # perform additional checks for parameters that require a specific structure (e.g lists of a given type)
    if name == "pos":
        for index, coord in enumerate(value):
            if not isinstance(coord, float):
                try:
                    value[index] = float(value[index])
                except:
                    raise ValueError(f"Invalid [{name}] input. Expected list[float].")

    if name in ["time", "sources"]:
        for element in value:
            if not isinstance(element, int):
                try:
                    element = int(element)
                except:
                    raise ValueError(f"Invalid [{name}] input. Expected list[int].")

    if name in ["colours", "plot_bands"]:
        for index, element in enumerate(value):
            if not isinstance(element, str):
                try:
                    value[index] = str(value[index])
                except:
                    raise ValueError(f"Invalid [{name}] input. Expected list[str].")

    if name == "overlays":
        supported_overlays = OverlayStruct().supportedOverlays

        if isinstance(value, list):
            for i, val in enumerate(value):
                value[i] = str(val).lower()
                if value[i] not in supported_overlays:
                    raise ValueError(f"Unsupported overlay '{val}'.")
        elif isinstance(value, dict):
            for key, val in value.items():
                if key not in supported_overlays:
                    raise ValueError(f"Unsupported overlay: '{val}'.")

    if name == "selection":
        supported_surveys = surveyInfo.defaultDataSurveys

        for entry in value:
            if not isinstance(entry, dict):
                raise ValueError("Invalid datatable selection. Expected list<dict>.")

        for element in value:
            if "kind" not in element:
                raise ValueError("Datatable selection missing key [kind].")

            if not isinstance(element["kind"], str):
                try:
                    element["kind"] = str(element["kind"]).lower()
                except:
                    raise ValueError(f"Invalid selection kind. Expected str, got {type(element['kind'])}.")

            if "survey" in element:
                if element["survey"] and not isinstance(element["survey"], str):
                    try:
                        element["survey"] = str(element["survey"])
                    except:
                        raise ValueError(f"Invalid selection survey. Expected str, got {type(element['survey'])}.")

            if element["kind"] == "atk_defaults":
                required_keys = ["kind", "surveys"]
                keys_to_type_check = ["surveys"]
            elif element["kind"] == "external":
                required_keys = ["kind", "survey", "parameters", "values", "errors", "notes"]
                keys_to_type_check = ["parameters", "values", "errors", "notes"]
            else:
                required_keys = ["kind", "survey", "parameters", "errors", "notes"]
                keys_to_type_check = ["parameters", "errors", "notes"]

            if not all(key in element for key in required_keys):
                raise ValueError(f"Datatable selection missing keys. Required keys: {required_keys}.")

            accepted_kinds = ["vizier", "atk_defaults", "external"]
            if element["kind"] not in accepted_kinds:
                raise ValueError(f"Invalid datatable selection kind. Accepted kinds: {accepted_kinds}.")

            for key in keys_to_type_check:
                if not isinstance(element[key], list):
                    try:
                        element[key] = [element[key]]
                    except:
                        raise ValueError(
                            f"Invalid value for key [{key}] in datatable selection. Expected list, got {type(element[key])}."
                        )

            if element["kind"] in ["atk_defaults"]:
                for survey in element["surveys"]:
                    if survey not in supported_surveys:
                        raise ValueError(f"Invalid selection survey. Supported surveys: {supported_surveys}.")

    if name == "dimensions":
        if "width" not in value or "height" not in value:
            raise ValueError("Grid dimensions missing required keys.")
        for key, val in value.items():
            if not isinstance(val, int):
                try:
                    value[key] = int(val)
                except:
                    raise ValueError("Grid dimensions must be integers.")

    if name == "panels":
        for plot_info in value:
            if (
                "name" not in plot_info
                or "width" not in plot_info
                or "height" not in plot_info
                or "figure" not in plot_info
            ):
                raise ValueError("datapage panels missing required keys.")
        if not isinstance(plot_info["name"], str):
            try:
                plot_info["name"] = str(plot_info["name"])
            except:
                raise ValueError(f"Plot names must be string, got {type(plot_info['name'])}.")
        if not isinstance(plot_info["width"], int):
            try:
                plot_info["width"] = int(plot_info["width"])
            except:
                raise ValueError(f"Plot widths must int, got {type(plot_info['width'])}.")
        if not isinstance(plot_info["height"], int):
            try:
                plot_info["height"] = int(plot_info["height"])
            except:
                raise ValueError(f"Plot heights must be int, got {type(plot_info['height'])}.")

    return value


def targeting_check(input):
    if input["pos"][0] and input["source"][0]:
        raise ValueError("Simultaneous source and pos input detected.")
    elif not input["pos"] and not input["source"]:
        raise ValueError("Source or pos input required.")


def check_inputs(inputs, label, check_targeting=False):
    """
    Checks inputs, using a list of lists in form:
    {name:[value,type],name:[value,type],...}
    """

    if not inputs:
        return inputs

    """ 
    Primary checks
    """
    corrected_inputs = {}
    for name, entry in inputs.items():
        value, target_type = entry
        if value is not None and value != "config":
            # automatic conversion to lists
            if name in ["sources", "columns", "colours", "plot_bands"] and not isinstance(value, list):
                value = [value]
            if name in ["overlays"] and not isinstance(value, (list, dict)):
                value = [value]
            corrected_value = check_type(name, value, target_type)
        else:
            corrected_value = value
        corrected_inputs[name] = corrected_value

    # check pos/source input configuration
    if check_targeting:
        targeting_check(inputs)

    """ 
    Additional (per-function) checks
    """
    if label == "query":
        kind = corrected_inputs["kind"]
        if kind not in toolInfo.supported_query_kinds:
            raise ValueError(f"Unsupported query kind '{kind}'. Accepted kinds: {toolInfo.supported_query_kinds}.")

        # check supported surveys
        survey = corrected_inputs["survey"]
        if corrected_inputs["kind"] not in ["data", "hrd", "bulkdata", "sed"]:
            if kind == "lightcurve" and survey == "gaia":
                survey = "gaia_lc"

            supported_surveys = getattr(surveyInfo, f"default{corrected_inputs['kind'].capitalize()}Surveys")
            if survey not in supported_surveys:
                raise ValueError(f"Unsupported survey in {kind} query. Supported surveys are: {supported_surveys}")

            if kind == "lightcurve" and survey == "gaia_lc":
                corrected_inputs["survey"] = "gaia"

        if kind == "data" and survey == "gaia_lc" and corrected_inputs["level"] != "internal":
            raise ValueError(f"Unsupported survey in {kind} query. Use kind=lightcurve for gaia lightcurve queries.")

        if kind == "image":
            if survey == "panstarrs":
                if corrected_inputs["size"] > 1500:
                    raise ValueError('Size too large. Maximum supported by panstarrs is 1500".')
                if not re.match("^[grizy]+$", corrected_inputs["band"]):
                    raise ValueError("Invalid panstarrs bands. Supported bands are [g,r,i,z,y].")
                corrected_inputs["size"] = int(corrected_inputs["size"])

            elif survey == "skymapper":
                if corrected_inputs["size"] > 600:
                    raise ValueError('Size too large. Maximum supported by skymapper is 600".')
                if not re.match("^[grizuv]+$", corrected_inputs["band"]):
                    raise ValueError("Invalid skymapper bands. Supported bands are [g,r,i,z,u,v].")
                else:
                    band = [corrected_inputs["band"]]
                    temp_string = ""
                    for i in range(0, len(band)):
                        temp_string += band[i] + ","
                    band = temp_string[:-1]

            elif survey == "dss":
                if corrected_inputs["size"] > 7200:
                    raise ValueError('Size too large. Maximum supported by dss is 7200".')
                if corrected_inputs["band"] != "g":
                    raise ValueError("Invalid dss band. Only g band is supported by DSS.")

    elif label == "search":
        if corrected_inputs["kind"] not in ["vizier", "simbad"]:
            raise ValueError(f"Invalid search kind '{corrected_inputs['kind']}'.")

    elif label == "editconfig":
        from ..Configuration.baseconfig import ConfigStruct

        config = ConfigStruct()
        config.read_config()
        supported_keys = config.supported_keys

        if corrected_inputs["key"] not in supported_keys:
            raise ValueError(f"Invalid config key '{corrected_inputs['key']}'")

    elif label == "plot":
        if corrected_inputs["data_kind"] == "lightcurve":
            plot_kind = corrected_inputs["kind"]
            if plot_kind and plot_kind not in ["lightcurve", "powspec", "phasefold"]:
                raise ValueError("Invalid light curve plot kind.")

            if plot_kind == "lightcurve":
                from ..Plotting.lightcurveplotting import SupportedColours

                colours = corrected_inputs["colours"]
                supported_colours = SupportedColours("green").supported_colours
                if colours:
                    for colour in colours:
                        if colour not in supported_colours:
                            raise ValueError(
                                f"Invalid light curve colour {colour}. Supported colours are: {supported_colours}"
                            )

                timeformat = corrected_inputs["timeformat"]
                if timeformat and timeformat not in ["reduced", "original"]:
                    raise ValueError("Invalid light curve time format.")

            elif plot_kind in ["powspec", "phasefold"]:
                method = corrected_inputs["method"]
                if method and method not in ["ls"]:
                    raise ValueError("Invalid time series method, accepted methods are: ls")

    elif label == "bin":
        binsize = corrected_inputs["binsize"]
        if binsize and binsize[-1] not in ["d", "h", "m"]:
            raise ValueError("Invalid binsize unit. Accepted units are: d (days), h (hours) and m (mins).")

    elif label == "crop":
        start, stop = corrected_inputs["start"], corrected_inputs["stop"]
        if start[-1] == "%" and not 0 <= float(start.rstrip("%")) <= 100:
            raise ValueError("Invalid crop 'start' percentage.")
        if stop[-1] == "%" and not 0 <= float(stop.rstrip("%")) <= 100:
            raise ValueError("Invalid crop 'stop' percentage.")

        if not (corrected_inputs["start"] and corrected_inputs["stop"]):
            raise ValueError("Invalid crop parameters. Requires 'start' and 'stop' as absolute values or percentages.")

    elif label == "correctpm":
        pos = corrected_inputs["pos"]
        source = corrected_inputs["source"]
        pmra, pmdec = corrected_inputs["pmra"], corrected_inputs["pmdec"]
        input_time, target_time = (corrected_inputs["input_time"], corrected_inputs["target_time"])
        input_survey, target_survey = (corrected_inputs["input_survey"], corrected_inputs["target_survey"])

        if pos and pmra and pmdec and input_time and target_time:
            pass
        elif pos and pmra and pmdec and input_survey and target_survey:
            pass
        elif source and input_time and target_time:
            pass
        elif source and target_time:
            pass
        elif source and target_survey:
            pass
        else:
            raise ValueError("Invalid correctpm input combination.")

    return list(corrected_inputs.values())
