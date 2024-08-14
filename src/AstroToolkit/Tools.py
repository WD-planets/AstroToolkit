"""
Contains the main tools for data fetching, plotting and analysis in ATK.
"""

from .Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()
newline = "\n"


# 'pos' has epoch of 2000 if given as input, 2016 if found as a result of source query. 'identifier' is therefore always J2000


def query(
    kind,
    survey=None,
    radius="config",
    pos=None,
    source=None,
    size="config",
    band="config",
    username="config",
    password="config",
    overlays="config",
    sources=None,
    level="external",
    raw=False,
):
    """query(kind,source/pos, *)

    Returns data of a given type from a given survey.

    :param str kind: Type of query to perform, from:

        - data: return survey data as listed in Vizier
        - phot: return only photometry from supported surveys
        - bulkphot: return photometry from all supported surveys
        - reddening: return reddening from a supported survey
        - lightcurve: return light curve data from a supported survey
        - image: return image data from a supported survey
        - hrd: return Gaia Hertzsprung-Russell diagram data
        - sed: return spectral energy distribution data from all supported surveys
        - spectrum: return spectrum data from a supported survey

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    Requires additional parameters depending on query type, and returns different data structures in each case.

    - kind = data

    :param str survey: Target survey, from:

        - any Vizier survey ID
        - gaia
        - panstarrs
        - skymapper
        - galex
        - rosat
        - sdss
        - wise
        - twomass (2MASS)
        - erosita

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`DataStruct`
    :rtype: class

    |

    - kind = phot

    :param str survey: Target survey, from:

        - gaia
        - panstarrs
        - skymapper
        - galex
        - sdss
        - wise
        - twomass (2MASS)

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`DataStruct`
    :rtype: class

    |

    - kind = bulkphot

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`DataStruct`
    :rtype: class

    |

    - kind = reddening

    :param str survey: Target survey, from:

        - stilism, requires source input, doesn't accept a radius
        - gdre, accepts source or pos input

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`DataStruct`
    :rtype: class

    |

    - kind = lightcurve

    :param str survey: Target survey, from:

        - atlas, requires username and password
        - ztf
        - crts
        - asassn
        - gaia
        - tess

    :param float, optional radius: Search radius in arcseconds, default: config
    :param bool, optional raw: Return raw data with no filtering, default: False
    :param str, optional username: ATLAS username, only required in ATLAS queries
    :param str, optional password: ATLAS password, only required in ATLAS queries

    :return: :ref:`LightcurveStruct`
    :rtype: class

    |

    - kind = image

    :param str survey: Target survey, from:

        - panstarrs, accepted bands = grizy
        - skymapper, accepted bands = grizuv
        - dss, accepted bands = g

    :param float, optional size: Size of image in arcseconds, default: config
    :param str, optional, band: All required image bands (as listed above) as a single string. E.g. 'grizy' for all panstarrs bands. Default = g
    :param list<str>, optional overlays: Required detection overlays from:

        - gaia
        - galex
        - wise
        - sdss
        - twomass  (2MASS)
        - skymapper
        - panstarrs
        - rosat
        - erosita
        - atlas
        - gaia_lc (Gaia light curve)
        - asassn
        - crts
        - ztf

    The bands used in these detection overlays are taken from the config. If you wish to use multiple bands from a single survey, detections can instead be requested using a dict. E.g. to use all Gaia magnitudes:

    .. code-block::  python

        overlays={'gaia':['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']}

    :return: :ref:`ImageStruct`
    :rtype: class

    |

    - kind = hrd

    :param list<int> sources: Target Gaia DR3 source(s)

    :return: :ref:`HrdStruct`
    :rtype: class

    |

    - kind = sed

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`SedStruct`
    :rtype: class

    |

    - kind = spectrum

    :param str survey: Target survey, from:

        - sdss

    :param float, optional radius: Search radius in arcseconds, default: config

    :return: :ref:`SpectrumStruct`
    :rtype: class

    |

    """

    from .Data.dataquery import query as data_query

    config.read_config()

    from .Input.input_validation import check_inputs

    if not isinstance(overlays, (list, dict)) and overlays != "config":
        overlays = [overlays]

    inputs = check_inputs(
        {
            "survey": survey,
            "radius": radius,
            "pos": pos,
            "source": source,
            "size": size,
            "band": band,
            "username": username,
            "password": password,
            "overlays": overlays,
            "sources": sources,
        },
        kind,
    )

    survey, radius, pos, source, size, band, username, password, overlays, sources = (
        inputs["survey"],
        inputs["radius"],
        inputs["pos"],
        inputs["source"],
        inputs["size"],
        inputs["band"],
        inputs["username"],
        inputs["password"],
        inputs["overlays"],
        inputs["sources"],
    )

    if radius == "config" and kind != "image" and kind != "hrd":
        radius = float(getattr(config, f"query_{kind}_radius"))
    if size == "config":
        size = int(config.query_image_size)
    if band == "config":
        band = config.query_image_band
    if overlays == "config":
        overlays = [str(config.query_image_overlays)]
    if username == "config":
        username = config.query_lightcurve_atlas_username
    if password == "config":
        password = config.query_lightcurve_atlas_password
    if kind == "hrd":
        survey = "Gaia"
        radius = None

    if (
        config.enable_notifications == "True"
        and level != "internal"
        and kind not in ["image", "bulkphot", "sed", "hrd"]
    ):
        print(
            f"Running {survey} {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}"
        )
    elif (
        config.enable_notifications == "True"
        and level != "internal"
        and kind == "image"
    ):
        print(
            f"Running {survey} {kind} query{newline}source = {source}{newline}pos = {pos}{newline}size = {size}{newline}"
        )
    elif (
        config.enable_notifications == "True"
        and level != "internal"
        and kind in ["bulkphot", "sed"]
    ):
        print(
            f"Running {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}"
        )
    elif (
        config.enable_notifications == "True"
        and level != "internal"
        and kind in ["hrd"]
    ):
        print(f"Running {kind} query{newline}sources = {sources}{newline}")

    if kind == "data":
        data = data_query(survey=survey, radius=radius, pos=pos, source=source)
        if source and survey == "gaia" and data.data:
            data.pos = [data.data["ra"][0], data.data["dec"][0]]
        elif source and survey != "gaia":
            gaia_data = data_query(survey="gaia", radius=radius, pos=pos, source=source)
            if gaia_data.data:
                data.pos = [gaia_data.data["ra"][0], gaia_data.data["dec"][0]]
            else:
                return data

    elif kind == "spectrum":
        from .Data.spectrumquery import query as spectrum_query

        data = spectrum_query(survey=survey, radius=radius, pos=pos, source=source)

    elif kind == "image":
        from .Data.imagequery import query as image_query

        if survey == "any":
            for survey in ["panstarrs", "skymapper", "dss"]:
                data = image_query(
                    survey=survey,
                    size=size,
                    band=band,
                    pos=pos,
                    source=source,
                    overlays=overlays,
                )
                if data.data:
                    break
        else:
            data = image_query(
                survey=survey,
                size=size,
                band=band,
                pos=pos,
                source=source,
                overlays=overlays,
            )

    elif kind == "lightcurve":
        from .Data.lightcurvequery import query as lightcurve_query

        data = lightcurve_query(
            survey=survey,
            radius=radius,
            pos=pos,
            source=source,
            username=username,
            password=password,
            raw=raw,
        )

    elif kind == "phot":
        from .Data.photquery import query as phot_query

        data = phot_query(pos=pos, source=source, radius=radius, survey=survey)
    elif kind == "bulkphot":
        from .Data.photquery import bulkphot_query

        data = bulkphot_query(pos=pos, source=source, radius=radius)
    elif kind == "sed":
        from .Data.sedquery import query as sed_query

        data = sed_query(pos=pos, source=source, radius=radius)
    elif kind == "reddening":
        from .Data.reddeningquery import query as reddening_query

        data = reddening_query(survey=survey, source=source, pos=pos, radius=radius)
    elif kind == "hrd":
        from .Data.hrdquery import gather_data

        data = gather_data(sources)
    else:
        raise Exception("Invalid kind passed to query.")

    from .Misc.identifier_generation import identifier_from_pos

    if hasattr(data, "source"):
        if data.source and data.data:
            gaia_data = data_query(survey="gaia", radius=radius, pos=pos, source=source)
            if gaia_data.data:
                ra, dec = gaia_data.data["ra2000"][0], gaia_data.data["dec2000"][0]
                data.identifier = identifier_from_pos([ra, dec])
        else:
            data.identifier = identifier_from_pos(data.pos)
    elif hasattr(data, "sources") and data.data:
        identifiers = []
        for source in data.sources:
            gaia_data = data_query(
                survey="gaia", radius=radius, pos=pos, source=source
            ).data
            ra, dec = gaia_data["ra2000"][0], gaia_data["dec2000"][0]
            identifiers.append(identifier_from_pos([ra, dec]))
        data.identifiers = identifiers
    else:
        data.identifier = identifier_from_pos(data.pos)

    from .FileHandling.file_naming import name_file

    fname = name_file(data)
    data.dataname = fname

    if kind not in ["data", "phot", "bulkphot", "reddening"]:
        from .FileHandling.file_naming import generate_plotname

        generate_plotname(data)

    return data


def correctpm(
    input_survey=None,
    target_survey=None,
    pos=None,
    source=None,
    input_time=None,
    target_time=None,
    pmra=None,
    pmdec=None,
):
    """correctpm(source/pos, *)
    Corrects a system's coordinates for proper motion between times or supported surveys.

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    Also additional parameters depending on whether a source or pos is used:

    - for source input, requires one of:

    :param str target_survey: any supported survey in any supported query type
    :param list<int> target_time: time to correct coordinates to in format [year,month]

    - for pos input, requires:

    :param list<int> input_time: epoch of supplied coordinates in format [year,month]
    :param list<int> target_time: time to correct coordinates to in format [year,month]
    :param float pmra: proper motion in right ascension in mas/yr
    :param float pmdec: proper motion in declination in mas/yr

    or

    :param str input_survey: any supported survey in any supported query type
    :param str target_survey: any supported survey in any supported query type


    :return: [right ascension, declination] in degrees
    :rtype: list<int>

    |

    """
    from .Input.input_validation import check_inputs

    inputs = check_inputs(
        {
            "input_survey": input_survey,
            "target_survey": target_survey,
            "pos": pos,
            "source": source,
            "input_time": input_time,
            "target_time": target_time,
            "pmra": pmra,
            "pmdec": pmdec,
        },
        "correctpm",
    )

    input_survey, target_survey, pos, source, input_time, target_time, pmra, pmdec = (
        inputs["input_survey"],
        inputs["target_survey"],
        inputs["pos"],
        inputs["source"],
        inputs["input_time"],
        inputs["target_time"],
        inputs["pmra"],
        inputs["pmdec"],
    )

    if source and target_time:
        from .Misc.pmcorrection import autocorrect_source

        corrected_pos = autocorrect_source(source=source, target_time=target_time)
    elif source and target_survey:
        from .Misc.pmcorrection import autocorrect_source

        corrected_pos = autocorrect_source(source=source, target_survey=target_survey)

    elif pos and input_time and target_time:
        from .Misc.pmcorrection import correctpm

        corrected_pos = correctpm(input_time, target_time, pos[0], pos[1], pmra, pmdec)
    elif pos and input_survey and target_survey:
        from .Misc.pmcorrection import autocorrect_survey

        corrected_pos = autocorrect_survey(
            input_survey, target_survey, ra=pos[0], dec=pos[1], pmra=pmra, pmdec=pmdec
        )
    else:
        raise Exception("Invalid input combination passed to correctpm.")

    return corrected_pos


def readdata(fname):
    """
    Reads data from a local file created by ATK, recreating the inital data structure.

    :param str fname: name of file from which to read

    :return: ATK Data Structure of same type as original data
    :rtype: class

    |

    """
    from .FileHandling.file_reading import read_local_file
    from .Input.input_validation import check_inputs

    inputs = check_inputs({"fname": fname}, "readdata")
    fname = inputs["fname"]

    config.read_config()
    if config.enable_notifications == "True":
        print(f"Recreating data from local storage: {fname}{newline}")

    from .FileHandling.file_naming import name_file

    struct = read_local_file(fname)
    struct.dataname = name_file(struct)

    if struct.kind not in ["data", "phot", "bulkphot", "reddening"]:
        from .FileHandling.file_naming import generate_plotname

        generate_plotname(struct)

    def savedata(self, fname=None):
        if fname:
            if self.kind != "image":
                if not fname.endswith(".csv"):
                    fname += ".csv"
            else:
                if not fname.endswith(".fits"):
                    fname += ".fits"
        else:
            fname = self.dataname

        from .FileHandling.file_writing import generate_local_file

        success = generate_local_file(self, fname)

        if success:
            config.read_config()
            if config.enable_notifications == "True":
                print(f"Saving data to local storage: {fname}{newline}")

        return fname

    import types

    struct.savedata = types.MethodType(savedata, struct)

    from .Data.data_printing import print_data

    struct.showdata = types.MethodType(print_data, struct)

    return struct


def search(kind, radius="config", source=None, pos=None):
    """search(kind,source/pos, *)

    Searches for a given target in Vizier or SIMBAD.

    :param str kind: where to search for target, from:

        - vizier
        - simbad

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees
    :param float, optional radius: radius of search in arcseconds, default = config

    :return: None

    |

    """

    from .Input.input_validation import check_inputs
    from .Misc.search import do_search

    inputs = check_inputs(
        {"kind": kind, "radius": radius, "source": source, "pos": pos}, "search"
    )
    kind, radius, source, pos = (
        inputs["kind"],
        inputs["radius"],
        inputs["source"],
        inputs["pos"],
    )

    config.read_config()
    if radius == "config":
        radius = float(config.search_radius)

    if config.enable_notifications == "True":
        print(
            f"Running {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}"
        )

    do_search(kind=kind, radius=radius, source=source, pos=pos)

    return None


def readfits(fname, columns):
    """
    Reads columns from a .fits file.

    :param str fname: name of file from which to read
    :param list<str>/str columns: name(s) of column(s) to read

    :return: Returned column data. E.g. if [ra, dec] requested, returns [[ra], [dec]]
    :rtype: list<list>

    |

    """
    from .Input.input_validation import check_inputs
    from .Misc.fitsfiles import get_columns

    inputs = check_inputs({"fname": fname, "columns": columns}, "readfits")
    fname, columns = inputs["fname"], inputs["columns"]

    config.read_config()
    if config.enable_notifications:
        print(f"Reading local .fits file: {fname}")

    return get_columns(filename=fname, columns=columns)


def deg2hms(pos):
    """
    Converts coordinates in degrees to HMS±DMS format.

    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: coordinates in HMS±DMS format
    :rtype: str

    |

    """
    from .Input.input_validation import check_inputs
    from .Misc.identifier_generation import identifier_from_pos

    inputs = check_inputs({"pos": pos}, "deg2hms")
    pos = inputs["pos"]

    return identifier_from_pos(pos, kind="conversion")


def hms2deg(pos):
    """
    Converts coordinates in HMS±DMS format to degrees.

    :param str pos: position of target in HMS±DMS format, i.e. HHMMSS.SS... ± DDMMSS.SS...

    :return: [right ascension, declination] in degrees
    :rtype: list<float>

    |

    """
    from .Input.input_validation import check_inputs
    from .Misc.coordinate_conversions import conv_hms_to_deg

    inputs = check_inputs({"pos": pos}, "hms2deg")
    pos = inputs["pos"]

    return conv_hms_to_deg(pos)


def tsanalysis(data):
    """
    Opens GUI for PyAOV time series analysis.

    :param class LightcurveStruct: Light curve data in format of ATK LightcurveStruct

    :return: None

    |

    """

    import os
    from pathlib import Path

    from .Timeseries.pyaov import pyaov
    from .Timeseries.pyaov.pyaov_interface import get_analysis

    # path = Path(pyaov.__file__).parent.absolute()
    # if str(path) not in os.environ["PATH"]:
    #    os.environ["PATH"] += str(path)

    get_analysis(struct=data, gui=True)
