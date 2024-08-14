import math

import pandas as pd


def ReadLocalData(fname):
    data = pd.read_csv(fname)

    atk_subkind = data["atk_subkind"][0]
    identifier, source, pos, survey, catalogue = (
        data["atk_identifier"][0],
        data["atk_source"][0],
        [float(data["atk_pos_ra"][0]), float(data["atk_pos_dec"][0])],
        data["atk_survey"][0],
        data["atk_catalogue"][0],
    )
    if math.isnan(float(source)):
        source = None
    if math.isnan(pos[0]) and math.isnan(pos[1]):
        pos = None
    if not isinstance(survey, str):
        survey = None
    if not isinstance(catalogue, str):
        catalogue = None

    data.drop(
        [
            "atk_identifier",
            "atk_source",
            "atk_pos_ra",
            "atk_pos_dec",
            "atk_survey",
            "atk_catalogue",
        ],
        inplace=True,
        axis=1,
    )

    data = data.to_dict("list")

    from ..Data.dataquery import DataStruct

    recreated_data = DataStruct(
        survey=survey,
        catalogue=catalogue,
        source=str(source),
        pos=pos,
        identifier=identifier,
        data=data,
        sub_kind=atk_subkind,
    )

    return recreated_data


def ReadLocalPhot(fname):
    recreated_data = ReadLocalData(fname)
    return recreated_data


def ReadLocalBulkphot(fname):
    import csv

    with open(fname, "r") as file:
        content = list(csv.reader(file, delimiter=","))

    surveys = []
    survey = []
    for _, line in enumerate(content):
        if line:
            survey.append(line)
        else:
            surveys.append(survey)
            survey = []

    for i, survey in enumerate(surveys):
        for j, _ in enumerate(survey):
            conv = lambda i: i or None
            surveys[i][j] = [conv(i) for i in surveys[i][j]]

    data_arr = []
    for survey in surveys:
        headers = survey.pop(0)
        df = pd.DataFrame(survey, columns=headers)
        data_arr.append(df)

    source, pos, identifier = (
        data_arr[0]["atk_source"][0],
        [data_arr[0]["atk_pos_ra"][0], data_arr[0]["atk_pos_dec"][0]],
        data_arr[0]["atk_identifier"][0],
    )

    if pos == [None, None]:
        pos = None

    from ..Data.dataquery import SurveyInfo

    supported_surveys = SurveyInfo().bulkphot_surveys

    recreated_data = {}
    for df in data_arr:
        survey = df["atk_survey"][0]
        df.drop(
            ["atk_source", "atk_pos_ra", "atk_pos_dec", "atk_survey", "atk_identifier"],
            inplace=True,
            axis=1,
        )
        data_dict = pd.DataFrame.to_dict(df, orient="list")
        recreated_data[survey] = data_dict

    for survey in supported_surveys:
        if survey not in recreated_data:
            recreated_data[survey] = None

    from ..Data.dataquery import DataStruct

    return DataStruct(
        survey="all",
        catalogue=None,
        source=str(source),
        pos=pos,
        identifier=identifier,
        data=recreated_data,
        sub_kind="bulkphot",
    )


def ReadLocalLightcurve(fname):
    data = pd.read_csv(fname, index_col=None, sep=",")
    bands = [x for x in dict.fromkeys(data["band"].tolist()) if x != "band"]
    data_arr = [data[data["band"] == band].reset_index(drop=True) for band in bands]

    recreated_data = []
    for element in data_arr:
        identifier, source, pos, survey, band = (
            element["atk_identifier"][0],
            element["atk_source"][0],
            [float(element["atk_pos_ra"][0]), float(element["atk_pos_dec"][0])],
            element["atk_survey"][0],
            element["band"][0],
        )

        ra, dec, mag, mag_err = (
            [float(x) for x in element["ra"].tolist()],
            [float(x) for x in element["dec"].tolist()],
            [float(x) for x in element["mag"].tolist()],
            [float(x) for x in element["mag_err"].tolist()],
        )

        if "hjd" in element.columns.values:
            time_type = "hjd"
        elif "mjd" in element.columns.values:
            time_type = "mjd"
        else:
            raise Exception("No/Invalid time type found in data. Expected mjd/hjd")

        time, time_ori = (
            [float(x) for x in element[time_type].tolist()],
            [float(x) for x in element[f"{time_type}_ori"].tolist()],
        )

        if math.isnan(float(source)):
            source = None
        if math.isnan(pos[0]) and math.isnan(pos[1]):
            pos = None

        recreated_data.append(
            {
                "band": band,
                "ra": ra,
                "dec": dec,
                time_type: time,
                f"{time_type}_ori": time_ori,
                "mag": mag,
                "mag_err": mag_err,
            }
        )

    from ..Data.dataquery import SurveyInfo

    lightcurve_bands = SurveyInfo().lightcurve_bands
    supported_bands = lightcurve_bands[survey]

    recreated_bands = [x["band"] for x in recreated_data]
    for band in supported_bands:
        if band not in recreated_bands:
            recreated_data.append(
                {
                    "band": band,
                    "ra": None,
                    "dec": None,
                    time_type: None,
                    f"{time_type}_ori": None,
                    "mag": None,
                    "mag_err": None,
                }
            )

    from ..Data.lightcurvequery import LightcurveStruct

    struct = LightcurveStruct(
        survey=survey,
        source=str(source),
        pos=pos,
        identifier=identifier,
        data=recreated_data,
    )

    return struct


def ReadLocalImage(fname):
    from astropy.io import fits
    from astropy.wcs import WCS

    image_data = fits.open(fname)[0]
    data, header = image_data.data, image_data.header

    survey, source, pos, identifier, image_focus, size, image_time, wcs = (
        header["atk_survey"],
        header["atk_source"],
        [float(header["atk_pos_ra"]), float(header["atk_pos_dec"])],
        header["atk_identifier"],
        [float(header["atk_image_focus_ra"]), float(header["atk_image_focus_dec"])],
        header["atk_size"],
        [header["atk_time_year"], header["atk_time_month"]],
        WCS(header),
    )

    # breaks when the end of the overlay entries has been reached, at which point i=number of data points. Can't use lists as not supported by fits.
    i = 0
    while True:
        try:
            overlay_len_check = header[f"atk_overlay_survey_{i}"]
            i += 1
        except:
            break
    overlay_len = i

    overlay_data = []
    for i in range(0, overlay_len):
        overlay_data_point = {
            "overlay_type": header[f"atk_overlay_type_{i}"],
            "marker_type": header[f"atk_overlay_marker_type_{i}"],
            "corrected": header[f"atk_overlay_corrected_{i}"],
            "ra": header[f"atk_overlay_ra_{i}"],
            "dec": header[f"atk_overlay_dec_{i}"],
            "marker_size": header[f"atk_overlay_marker_size_{i}"],
            "colour": header[f"atk_overlay_colour_{i}"],
            "mag_name": header[f"atk_overlay_mag_name_{i}"],
            "survey": header[f"atk_overlay_survey_{i}"],
        }
        overlay_data.append(overlay_data_point)

    from ..Data.imagequery import ImageStruct

    recreated_data = ImageStruct(
        survey=survey,
        source=source,
        pos=pos,
        identifier=identifier,
        data={
            "image_data": data,
            "image_header": header,
            "size": size,
            "image_time": image_time,
            "wcs": wcs,
            "image_focus": image_focus,
            "overlay": overlay_data,
        },
    )
    return recreated_data


def ReadLocalReddening(fname):
    recreated_data = ReadLocalData(fname)
    return recreated_data


def ReadLocalSed(fname):
    data = pd.read_csv(fname, index_col=None, sep=",")
    surveys = [
        x for x in dict.fromkeys(data["sed_survey"].tolist()) if x != "sed_survey"
    ]
    data_arr = [
        data[data["sed_survey"] == survey].reset_index(drop=True) for survey in surveys
    ]

    recreated_data = []
    for element in data_arr:
        survey = element["sed_survey"][0]
        source = element["atk_source"][0]
        pos = [element["atk_pos_ra"][0], element["atk_pos_dec"][0]]
        identifier = element["atk_identifier"][0]

        data_dict = {
            "survey": survey,
            "wavelength": element["wavelength"].tolist(),
            "flux": element["flux"].tolist(),
            "flux_rel_err": element["flux_rel_err"].tolist(),
        }
        recreated_data.append(data_dict)

    from ..Data.sedquery import SedStruct

    return SedStruct(source=source, pos=pos, identifier=identifier, data=recreated_data)


def ReadLocalSpectrum(fname):
    data = pd.read_csv(fname, index_col=None, sep=",")

    pos_ra, pos_dec = data["atk_pos_ra"][0], data["atk_pos_dec"][0]
    recreated_data = {
        "wavelength": data["wavelength"].tolist(),
        "flux": data["flux"].tolist(),
    }

    from ..Data.spectrumquery import SpectrumStruct

    return SpectrumStruct(
        survey=data["atk_survey"][0],
        source=data["atk_source"][0],
        pos=[pos_ra, pos_dec],
        identifier=data["atk_identifier"][0],
        data=recreated_data,
    )


def ReadLocalHrd(fname):
    data = pd.read_csv(fname, index_col=None, sep=",")
    recreated_data = {"bp-rp": data["bp-rp"], "absg": data["absg"]}

    from ..Data.hrdquery import HrdStruct

    return HrdStruct(
        survey=data["atk_survey"],
        sources=data["atk_sources"],
        identifiers=data["atk_identifiers"],
        data=recreated_data,
    )


def read_local_file(fname):
    if fname.endswith(".csv"):
        data = pd.read_csv(fname)
        atk_kind = data["atk_kind"][0]
    elif fname.endswith(".fits"):
        atk_kind = "image"

    data = globals()[f"ReadLocal{atk_kind.capitalize()}"](fname)

    return data
