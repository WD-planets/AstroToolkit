"""
Contains the main tools for data fetching, plotting and analysis in ATK.
"""

from .Configuration.baseconfig import ConfigStruct
from .Configuration.epochs import EpochStruct
from .Data.dataquery import DataStruct
from .Data.hrdquery import HrdStruct
from .Data.imagequery import ImageStruct
from .Data.lightcurvequery import LightcurveStruct
from .Data.sedquery import SedStruct
from .Data.spectrumquery import SpectrumStruct
from .Input.input_validation import check_inputs

config = ConfigStruct()
config.read_config()
newline = "\n"

# 'pos' has epoch of 2000 if given as input, 2016 if found as a result of source query. 'identifier' is therefore always J2000

epochs = EpochStruct().get_epochs()


def query(
    kind: str,
    survey: str = None,
    radius="config",
    pos: list[float] = None,
    source: int = None,
    size: float = "config",
    band: str = "config",
    username: str = "config",
    password: str = "config",
    overlays: list[str] = "config",
    sources: list[int] = None,
    level: str = "external",
    raw: bool = False,
    check_exists: str = None,
    retry: int = 1,
) -> DataStruct | HrdStruct | LightcurveStruct | ImageStruct | SedStruct | SpectrumStruct:
    """query(kind, source/pos, check_exists, **kwargs)
    Returns a :ref:`data structure <Data Structures>` of a given type from a given survey. Accepted types are:

    - :ref:`data <data-query>`: returns catalogue data as listed in Vizier
    - :ref:`bulkdata <bulkdata-query>`: returns data from all supported catalogues (including user-defined :ref:`aliases <Adding Catalogue Aliases>`)
    - :ref:`reddening <reddening-query>`: returns reddening from a supported survey
    - :ref:`lightcurve <lightcurve-query>`: returns light curve data from a supported survey
    - :ref:`image <image-query>`: returns image data from a supported survey
    - :ref:`hrd <hrd-query>`: returns Gaia Hertzsprung-Russell diagram data
    - :ref:`sed <sed-query>`: returns spectral energy distribution data from all supported surveys
    - :ref:`spectrum <spectrum-query>`: returns spectrum data from a supported survey

    |

    The type of query to be performed is chosen via the 'kind' parameter:

    :param kind: Type of query to perform, as listed above
    :type kind: str

    :func:`query()` requires additional parameters depending on the kind of query being performed, and returns different data structures in each case.

    |

    .. _data-query:
    .. rubric:: Data Queries
        :heading-level: 3

    :param survey: Target survey, from :ref:`supported surveys <Data Surveys>`
    :type survey: str
    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`query_data_radius <cfg_query_data_radius>` config key
    :type radius: float, optional
    :param check_exists: Path to check for existing data. If a file is found, data is generated without having to run a query. If a file is not found, query will go ahead and the resulting data will be saved to the requested Path for future executions. Defaults to None (i.e. this functionality is disabled)
    :type check_exists: bool, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`

    |

    .. _bulkdata-query:
    .. rubric:: Bulkdata Queries
        :heading-level: 3

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`query_bulkdata_radius <cfg_query_bulkdata_radius>` config key
    :type radius: float, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`
    :rtype: class

    |

    .. _reddening-query:
    .. rubric:: Reddening Queries
        :heading-level: 3

    :param str survey: Target survey, from:

        - stilism - requires source input, **doesn't accept a radius**
        - gdre - accepts source or pos input

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`query_reddening_radius <cfg_query_reddening_radius>` config key
    :type radius: float, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional


    :return: :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`

    |

    .. _lightcurve-query:
    .. rubric:: Light Curve Queries
        :heading-level: 3

    :param survey: Target survey, from :ref:`supported surveys <Light Curve Surveys>`
    :type survey: str

    :param radius: Search radius in arcseconds, default given by :ref:`query_lightcurve_radius <cfg_query_lightcurve_radius>` config key
    :param raw: Return raw data with no filtering, defaults to False
    :type raw: bool, optional
    :param username: ATLAS username, only required in ATLAS queries. Default given by :ref:`query_lightcurve_atlas_username <cfg_query_lightcurve_atlas_username>` config key
    :type usename: str, optional
    :param password: ATLAS password, only required in ATLAS queries. Default given by :ref:`query_lightcurve_atlas_password <cfg_query_lightcurve_atlas_password>` config key
    :type password: str, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>`

    |

    .. _image-query:
    .. rubric:: Image Queries
        :heading-level: 3

    :param survey: Target survey, from :ref:`supported surveys <Image Surveys>`
    :type survey: str

    :param size: Size of image in arcseconds, default given by :ref:`query_image_size <cfg_query_image_size>` config key
    :param size: float, optional
    :param band: All required image bands as a single string (see :ref:`supported bands <Image Surveys>`). E.g. 'grizy' for all panstarrs bands. Defaults to 'g'
    :type band: str, optional
    :param overlays: Required detection overlays. Accepts any :ref:`data survey <Data Surveys>` or :ref:`light curve survey <Light Curve Surveys>`
    :type overlays: list<str>, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    **Note:** use 'gaia_lc' for Gaia light curve overlays, and 'gaia' for gaia detection overlays

    The default band in the detection overlay of a given survey is taken from the relevant :ref:`[survey]_overlay_mag <Image Overlay Settings>` config key. If you wish to use multiple bands from a single survey, detections can instead be requested using a dictionary.

    E.g. to request an overlay that includes all Gaia magnitudes:

    .. code-block::  python

        overlays={'gaia':['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']}

    :return: :class:`ImageStruct <AstroToolkit.Data.imagequery.ImageStruct>`

    |

    .. _hrd-query:
    .. rubric:: HRD Queries
        :heading-level: 3

    :param sources: List of target Gaia DR3 sources
    :type sources: list<int>
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`HrdStruct <AstroToolkit.Data.hrdquery.HrdStruct>`

    |

    .. _sed-query:
    .. rubric:: SED Queries
        :heading-level: 3

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`query_sed_radius <cfg_query_sed_radius>` config key
    :type radius: float, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`SedStruct <AstroToolkit.Data.sedquery.SedStruct>`

    |

    .. _spectrum-query:
    .. rubric:: Spectrum Queries
        :heading-level: 3

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param survey: Target survey, from :ref:`supported surveys <Spectrum Surveys>`
    :type survey: str
    :param radius: Search radius in arcseconds, default given by :ref:`query_spectrum_radius <cfg_query_spectrum_radius>` config key
    :type radius: float, optional
    :param retry: Re-attempt query if no data is returned. This can help to ensure that no data is missed due to e.g. brief server interruptions when performing queries for a large number of targets. Defaults to 1 (i.e. only one attempt is made).
    :type retry: int, optional

    :return: :class:`SpectrumStruct <AstroToolkit.Data.spectrumquery.SpectrumStruct>`

    |

    """
    retry_count, data_found = 0, False

    save_data = False
    if check_exists:
        import os
        from pathlib import Path

        try:
            path = Path(check_exists)
        except:
            raise ValueError("Invalid path.")

        base_file = os.path.basename(path)
        if not base_file.endswith(".fits"):
            base_file += ".fits"
        path = Path(os.path.join(path.parent.absolute(), base_file))

        if path.is_file():
            print(f"Note: existing {kind} data found in {path}")
            return readdata(path)
        else:
            save_data = True

    from .Data.dataquery import query as data_query

    config.read_config()

    # gets config values (if required by given query kind)
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

    # check inputs
    corrected_inputs = check_inputs(
        {
            "kind": [kind, str],
            "survey": [survey, str],
            "radius": [radius, float],
            "pos": [pos, list],
            "source": [source, int],
            "size": [size, float],
            "band": [band, str],
            "username": [username, str],
            "password": [password, str],
            "overlays": [overlays, (list, dict)],
            "sources": [sources, list],
            "level": [level, str],
            "raw": [raw, bool],
        },
        label="query",
        check_targeting=True,
    )

    # get validated inputs
    (kind, survey, radius, pos, source, size, band, username, password, overlays, sources, level, raw) = (
        corrected_inputs
    )

    # Handles notification output depending on query kind
    if config.enable_notifications and level != "internal":
        if kind == "image":
            print(
                f"{newline}Running {survey} {kind} query{newline}source = {source}{newline}pos = {pos}{newline}size = {size}{newline}"
            )
        elif kind in ["bulkdata", "sed"]:
            print(
                f"{newline}Running {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}"
            )
        elif kind in ["hrd"]:
            print(f"{newline}Running {kind} query{newline}sources = {sources}{newline}")
        else:
            print(
                f"{newline}Running {survey} {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}"
            )

    # perform queries based on kind
    while retry_count < retry and not data_found:
        if retry_count > 0:
            print(f"{newline}Retrying... (attempt #{retry_count + 1}){newline}")

        if kind == "data":
            data = data_query(survey=survey, radius=radius, pos=pos, source=source)

        elif kind == "spectrum":
            from .Data.spectrumquery import query as spectrum_query

            data = spectrum_query(survey=survey, radius=radius, pos=pos, source=source)

        elif kind == "image":
            from .Data.imagequery import query as image_query

            data = image_query(survey=survey, size=size, band=band, pos=pos, source=source, overlays=overlays)

        elif kind == "lightcurve":
            from .Data.lightcurvequery import query as lightcurve_query

            data = lightcurve_query(
                survey=survey, radius=radius, pos=pos, source=source, username=username, password=password, raw=raw
            )

        elif kind == "bulkdata":
            from .Data.bulkquery import bulkdata_query

            data = bulkdata_query(pos=pos, source=source, radius=radius)

        elif kind == "sed":
            from .Data.sedquery import query as sed_query

            data = sed_query(pos=pos, source=source, radius=radius)

        elif kind == "reddening":
            from .Data.reddeningquery import query as reddening_query

            data = reddening_query(survey=survey, source=source, pos=pos, radius=radius)

        elif kind == "hrd":
            from .Data.hrdquery import gather_data

            data = gather_data(sources)

        if kind in ["lightcurve"]:
            for band in data.data:
                if band["mag"]:
                    data_found = True
            if not data_found:
                retry_count += 1
        else:
            if not data.data:
                retry_count += 1
            else:
                data_found = True

    from .Misc.identifier_generation import identifier_from_pos

    # generates identifiers in HRD queries
    if hasattr(data, "sources"):
        identifiers = []
        for position in data.positions:
            identifiers.append(identifier_from_pos(position))

        data.identifiers = identifiers
    else:
        data.identifier = identifier_from_pos(data.pos)

    from .FileHandling.file_naming import name_file

    # generates data names
    fname = name_file(data)
    data.dataname = fname

    # generates plot names
    if kind not in ["data", "bulkdata", "reddening"]:
        from .FileHandling.file_naming import generate_plotname

        generate_plotname(data)

    if save_data:
        data.savedata(check_exists)

    return data


def correctpm(
    input_survey: str = None,
    target_survey: str = None,
    pos: list[float] = None,
    source: int = None,
    input_time: list[int] = None,
    target_time: list[int] = None,
    pmra: float = None,
    pmdec: float = None,
    check_success: bool = False,
) -> list[float]:
    """correctpm(source/pos, **kwargs)
    Corrects coordinates for proper motion between times or supported surveys.

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>

    **Note:** Additional parameters are required depending on whether a source or pos is used:

    |

    .. rubric:: Source Input
        :heading-level: 3

    Requires one of:

    :param target_survey: any :ref:`supported survey <Supported Surveys>`
    :type target_survey: str
    :param target_time: target time in format [year,month]
    :type target_time: list<int>

    .. rubric:: Pos Input
        :heading-level: 3

    Requires:

    :param input_time: epoch of supplied coordinates in format [year,month]
    :type input_time: list<int>
    :param target_time: target time in format [year,month]
    :type target_time: list<int>
    :param pmra: proper motion in right ascension in mas/yr
    :type pmra: float
    :param pmdec: proper motion in declination in mas/yr
    :type pmdec: float

    or

    :param input_survey: any :ref:`supported survey <Supported Surveys>`
    :type input_survey: str
    :param target_survey: any :ref:`supported survey <Supported Surveys>`
    :type target_survey: str
    :param pmra: proper motion in right ascension in mas/yr
    :type pmra: float
    :param pmdec: proper motion in declination in mas/yr
    :type pmdec: float

    :return: [right ascension, declination] in degrees
    :rtype: list<int>

    |

    """

    corrected_inputs = check_inputs(
        {
            "input_survey": [input_survey, str],
            "target_survey": [target_survey, str],
            "pos": [pos, list],
            "source": [source, int],
            "input_time": [input_time, list],
            "target_time": [target_time, list],
            "pmra": [pmra, float],
            "pmdec": [pmdec, float],
            "check_success": [check_success, bool],
        },
        label="correctpm",
        check_targeting=True,
    )

    (input_survey, target_survey, pos, source, input_time, target_time, pmra, pmdec, check_success) = corrected_inputs

    if source and target_time:
        from .Misc.pmcorrection import autocorrect_source

        corrected_pos, success = autocorrect_source(source=source, target_time=target_time, check_success=True)
    elif source and target_survey:
        from .Misc.pmcorrection import autocorrect_source

        corrected_pos, success = autocorrect_source(source=source, target_survey=target_survey, check_success=True)

    elif pos and input_time and target_time:
        from .Misc.pmcorrection import correctpm

        corrected_pos, success = correctpm(input_time, target_time, pos[0], pos[1], pmra, pmdec, check_success=True)
    elif pos and input_survey and target_survey:
        from .Misc.pmcorrection import autocorrect_pos

        corrected_pos, success = autocorrect_pos(
            input_survey, target_survey, ra=pos[0], dec=pos[1], pmra=pmra, pmdec=pmdec, check_success=True
        )
    else:
        raise ValueError("Invalid input combination passed to correctpm.")

    if check_success:
        return corrected_pos, success
    else:
        return corrected_pos


def readdata(fname: str = None) -> DataStruct | HrdStruct | LightcurveStruct | ImageStruct | SedStruct | SpectrumStruct:
    """readdata(fname)
    Reads data from a local file created by ATK, recreating the original data structure. If no file name is provided, a file dialogue will open in which a file may be selected.

    :param fname: name of file from which to read
    :type fname: str, optional

    :return: :ref:`ATK Data Structure <Data Structures>`

    |

    """
    from .FileHandling.file_reading import read_local_file

    corrected_inputs = check_inputs({"fname": [fname, str]}, "readdata")
    fname = corrected_inputs[0]

    config.read_config()
    if config.enable_notifications:
        print(f"Recreating data from local storage: {fname}")

    from .FileHandling.file_naming import name_file

    if not fname:
        from Utility import openFileDialogue

        fname = openFileDialogue()

    struct = read_local_file(fname)
    struct.dataname = name_file(struct)

    if struct.kind not in ["data", "bulkdata", "reddening"]:
        from .FileHandling.file_naming import generate_plotname

        generate_plotname(struct)

    return struct


def search(kind: str, radius: float = "config", source: int = None, pos: float = None) -> None:
    """search(kind,source/pos, *)

    Searches for a given target in Vizier or SIMBAD.

    :param kind: where to search for target, from 'vizier', 'simbad'
    :type kind: str

    :param source: Target GAIA DR3 Source ID
    :type source: int
    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>
    :param radius: Search radius in arcseconds, default given by :ref:`search_radius <cfg_search_radius>` config key
    :type radius: float, optional

    :return: None

    |

    """

    from .Misc.search import do_search

    config.read_config()
    if radius == "config":
        radius = float(config.search_radius)

    corrected_inputs = check_inputs(
        {"kind": [kind, str], "radius": [radius, float], "source": [source, int], "pos": [pos, list]},
        "search",
        check_targeting=True,
    )
    kind, radius, source, pos = corrected_inputs

    if config.enable_notifications:
        print(f"Running {kind} query{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}")

    do_search(kind=kind, radius=radius, source=source, pos=pos)

    return None


def deg2hms(pos):
    """
    Converts coordinates in degrees to HMS±DMS format.

    :param pos: Position [right ascension, declination] in degrees
    :type pos: list<float>

    :return: coordinates in HMS±DMS format
    :rtype: str

    |

    """
    from .Misc.identifier_generation import identifier_from_pos

    corrected_inputs = check_inputs({"pos": [pos, list]}, "deg2hms")
    pos = corrected_inputs["pos"]

    return identifier_from_pos(pos, kind="conversion")


def hms2deg(identifier):
    """
    Converts coordinates in HMS±DMS format to degrees.

    :param str pos: position of target in HMS±DMS format, i.e. HHMMSS.SS...±DDMMSS.SS...

    :return: [right ascension, declination] in degrees
    :rtype: list<float>

    |

    """
    from .Misc.coordinate_conversions import conv_hms_to_deg

    corrected_inputs = check_inputs({"identifier": [identifier, str]}, "hms2deg")
    identifier = corrected_inputs["identifier"]

    return conv_hms_to_deg(identifier)


def readfits(fname, columns):
    """
    Reads columns from a .fits file.

    :param fname: name of file from which to read
    :type fname: str
    :param columns: names of column(s) to read
    :type columns: list<str>

    :return: Returned column data. E.g. if [ra, dec] requested, returns [[ra], [dec]]
    :rtype: list<list>

    |

    """
    from .Input.input_validation import check_inputs
    from .Misc.fitsfiles import get_columns

    corrected_inputs = check_inputs({"fname": [fname, str], "columns": [columns, list]}, "readfits")
    fname, columns = corrected_inputs

    config.read_config()
    if config.enable_notifications:
        print(f"Reading local .fits file: {fname}")

    return get_columns(filename=fname, columns=columns)


"""
def getspectype(sources):
    from .Misc.estimate_spectral_type import get_spectral_types

    corrected_inputs = check_inputs({"sources": [sources, list]}, "getspectype")
    sources = corrected_inputs[0]

    return get_spectral_types(sources)
"""
