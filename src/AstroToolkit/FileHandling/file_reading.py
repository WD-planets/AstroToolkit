from astropy.io import fits
from astropy.table import Table

from ..PackageInfo import SurveyInfo

surveyInfo = SurveyInfo()


def readHeader(hdul, *args):
    header = hdul[0].header
    values = []
    for entry in args:
        values.append(header[f"atk_{entry}"])
    return values


def checkTargeting(source, pos):
    if not source:
        source = None
    if not pos[0] and not pos[1]:
        pos = None
    return source, pos


def ReadLocalData(hdul):
    data = Table.read(hdul[1], format="fits").to_pandas()

    subkind, identifier, source, pos_ra, pos_dec, survey, catalogue, trace = readHeader(
        hdul, "subkind", "identifier", "source", "pos_ra", "pos_dec", "survey", "catalogue", "t"
    )
    pos = [pos_ra, pos_dec]
    source, pos = checkTargeting(source, pos)

    if not isinstance(survey, str):
        survey = None
    if not isinstance(catalogue, str):
        catalogue = None

    data = data.to_dict("list")

    from ..Data.dataquery import DataStruct

    recreated_data = DataStruct(
        survey=survey,
        catalogue=catalogue,
        source=source,
        pos=pos,
        identifier=identifier,
        data=data,
        sub_kind=subkind,
        trace=trace,
    )

    return recreated_data


def ReadLocalBulkdata(hdul):
    data = {}
    for hdu in hdul[1:]:
        survey = hdu.header["atk_survey"]
        data[survey] = Table.read(hdu, format="fits").to_pandas().to_dict("list")

    for survey in data:
        if not data[survey]:
            data[survey] = None

    from ..Data.dataquery import DataStruct

    source, pos_ra, pos_dec, identifier, trace = readHeader(hdul, "source", "pos_ra", "pos_dec", "identifier", "t")
    pos = [pos_ra, pos_dec]
    source, pos = checkTargeting(source, pos)

    recreated_data = DataStruct(
        survey="all",
        catalogue=None,
        source=source,
        pos=pos,
        identifier=identifier,
        data=data,
        sub_kind="bulkdata",
        trace=trace,
    )

    return recreated_data


def ReadLocalLightcurve(hdul):
    from ..Data.lightcurvequery import LightcurveStruct

    survey, source, pos_ra, pos_dec, identifier, trace = readHeader(
        hdul, "survey", "source", "pos_ra", "pos_dec", "identifier", "t"
    )
    pos = [pos_ra, pos_dec]
    source, pos = checkTargeting(source, pos)

    empty_band = {"ra": None, "dec": None, "mjd": None, "mag": None, "mag_err": None}

    data = []
    for hdu in hdul[1:]:
        entry = {}
        entry["band"] = hdu.header["atk_band"]
        band = Table.read(hdu, format="fits").to_pandas().to_dict("list")
        if band:
            for key, val in band.items():
                entry[key] = val
        else:
            for key, val in empty_band.items():
                entry[key] = val

        data.append(entry)

    recreated_data = LightcurveStruct(
        survey=survey, source=str(source), pos=pos, identifier=identifier, data=data, trace=trace
    )

    return recreated_data


def ReadLocalImage(hdul):
    from astropy.wcs import WCS

    image_data = hdul[0]
    data = image_data.data
    (
        survey,
        source,
        pos_ra,
        pos_dec,
        identifier,
        image_focus_ra,
        image_focus_dec,
        size,
        image_time_year,
        image_time_month,
        trace,
    ) = readHeader(
        hdul,
        "survey",
        "source",
        "pos_ra",
        "pos_dec",
        "identifier",
        "image_focus_ra",
        "image_focus_dec",
        "image_size",
        "image_time_year",
        "image_time_month",
        "t",
    )
    pos = [pos_ra, pos_dec]
    image_focus = [image_focus_ra, image_focus_dec]
    image_time = [image_time_year, image_time_month]

    wcs = WCS(image_data.header)

    overlay = []
    if len(hdul) > 1:
        for hdu in hdul[1:]:
            overlay_data = Table.read(hdu, format="fits").to_pandas()
            overlay_length = len(overlay_data)
            overlay_data = overlay_data.to_dict("list")
            overlay_keys = list(overlay_data.keys())

            survey_data = []
            for i in range(0, overlay_length):
                entry = {}
                for key in overlay_keys:
                    entry[key] = overlay_data[key][i]
                survey_data.append(entry)
            overlay += survey_data

    from ..Data.imagequery import ImageStruct

    recreated_data = ImageStruct(
        survey=survey,
        source=source,
        pos=pos,
        identifier=identifier,
        data={
            "image_data": data,
            "image_header": image_data.header,
            "size": size,
            "image_time": image_time,
            "wcs": wcs,
            "image_focus": image_focus,
            "overlay": overlay,
        },
        trace=trace,
    )

    return recreated_data


def ReadLocalReddening(fname):
    recreated_data = ReadLocalData(fname)
    return recreated_data


def ReadLocalSed(hdul):
    data = Table.read(hdul[1], format="fits").to_pandas()

    surveys = list(dict.fromkeys(data["survey"].tolist()))
    data_arr = [data[data["survey"] == survey].reset_index(drop=True) for survey in surveys]

    source, pos_ra, pos_dec, identifier, trace = readHeader(hdul, "source", "pos_ra", "pos_dec", "identifier", "t")
    pos = [pos_ra, pos_dec]

    recreated_data = []
    for element in data_arr:
        survey = element["survey"][0]

        data_dict = {
            "survey": survey,
            "mag_name": element["mag_name"].tolist(),
            "wavelength": element["wavelength"].tolist(),
            "flux": element["flux"].tolist(),
            "flux_rel_err": element["flux_rel_err"].tolist(),
        }
        recreated_data.append(data_dict)

    from ..Data.sedquery import SedStruct

    recreated_data = SedStruct(source=source, pos=pos, identifier=identifier, data=recreated_data, trace=trace)

    return recreated_data


def ReadLocalSpectrum(hdul):
    data = Table.read(hdul[1], format="fits").to_pandas()

    survey, source, pos_ra, pos_dec, identifier, trace = readHeader(
        hdul, "survey", "source", "pos_ra", "pos_dec", "identifier", "t"
    )

    recreated_data = {"wavelength": data["wavelength"].tolist(), "flux": data["flux"].tolist()}

    from ..Data.spectrumquery import SpectrumStruct

    return SpectrumStruct(
        survey=survey, source=source, pos=[pos_ra, pos_dec], identifier=identifier, data=recreated_data, trace=trace
    )


def ReadLocalHrd(hdul):
    data = Table.read(hdul[1], format="fits").to_pandas()

    survey, traces = readHeader(hdul, "survey", "t")
    sources = data["sources"].tolist()
    identifiers = data["identifiers"].tolist()
    positions_ra = data["position_ra"].tolist()
    positions_dec = data["position_dec"].tolist()
    positions = []
    for ra, dec in zip(positions_ra, positions_dec):
        positions.append([ra, dec])

    data.drop(columns=["sources", "identifiers"], inplace=True)

    recreated_data = {"bp-rp": data["bp-rp"].tolist(), "absg": data["absg"].tolist()}

    from ..Data.hrdquery import HrdStruct

    return HrdStruct(
        survey=survey, sources=sources, identifiers=identifiers, data=recreated_data, positions=positions, traces=traces
    )


def read_local_file(fname):
    hdul = fits.open(fname)
    header = hdul[0].header
    if "atk_kind" in header:
        atk_kind = header["atk_kind"]
    else:
        raise Exception("No atk_kind found in fits header.")

    data = globals()[f"ReadLocal{atk_kind.capitalize()}"](hdul)

    return data
