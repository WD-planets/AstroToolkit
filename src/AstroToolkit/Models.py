"""
This module provides empty ATK data structures, allowing ATK methods to be used on external data. These structures will initialise with all attributes set to None, unless a **source** or **position** is given in which case some attributes are automatically generated.

**Note:** the savedata method is not currently supported for custom data structures.
"""

from __future__ import annotations

import types

from .FileHandling.file_naming import generate_plotname, name_file
from .Input.input_validation import check_inputs
from .Misc.identifier_generation import identifier_from_pos
from .Tools import query


def savedata_disabled(_):
    raise Exception("savedata method is not supported for custom data models.")


def CustomDataStruct(source: int = None, pos: list[float] = None) -> CustomDataStruct:
    """CustomDataStruct(source/pos)
    Generates an empty :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`, as returned by :ref:`data queries <data-query>`.

    Optionally, the pos/source to which the data structure will refer may be passed. Some attributes will then be automatically filled. If a source or pos is not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: :class:`DataStruct <AstroToolkit.Data.dataquery.DataStruct>`
    :rtype: class

    |

    """

    from .Data.dataquery import DataStruct as BaseDataStruct

    corrected_inputs = check_inputs({"source": [source, int], "pos": [pos, list]}, "customdatastruct")
    source, pos = corrected_inputs

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia").data
        pos = [gaia_data["ra2000"][0], gaia_data["dec2000"][0]]
        identifier = identifier_from_pos(pos)
    elif pos:
        identifier = identifier_from_pos(pos)
    else:
        identifier = None

    CustomDataStruct = BaseDataStruct(
        survey=None, catalogue=None, source=source, pos=pos, identifier=identifier, data=None
    )

    if source or pos:
        fname = name_file(CustomDataStruct)
        CustomDataStruct.dataname = fname
    else:
        CustomDataStruct.dataname = None

    CustomDataStruct.subkind = None

    CustomDataStruct.savedata = types.MethodType(savedata_disabled, CustomDataStruct)

    return CustomDataStruct


def CustomLightcurveStruct(source: int = None, pos: list[float] = None) -> CustomLightcurveStruct:
    """CustomLightcurveStruct(source/pos)
    Generates an empty :Class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>`, as returned by :ref:`lightcurve queries <lightcurve-query>`.

    Optionally, the pos/source to which the data structure will refer may be passed. Some attributes will then be automatically filled. If a source or pos is not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: :class:`LightcurveStruct <AstroToolkit.Data.lightcurvequery.LightcurveStruct>`
    :rtype: class

    |

    """

    from .Data.lightcurvequery import LightcurveStruct as BaseLightcurveStruct

    corrected_inputs = check_inputs({"source": [source, int], "pos": [pos, list]}, "customlightcurvestruct")
    source, pos = corrected_inputs

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia").data
        pos = [gaia_data["ra2000"][0], gaia_data["dec2000"][0]]
        identifier = identifier_from_pos(pos)
    elif pos:
        identifier = identifier_from_pos(pos)
    else:
        identifier = None

    CustomLightcurveStruct = BaseLightcurveStruct(survey=None, source=source, pos=pos, identifier=identifier, data=None)

    if source or pos:
        fname = name_file(CustomLightcurveStruct)
        CustomLightcurveStruct.dataname = fname
        generate_plotname(CustomLightcurveStruct)
    else:
        CustomLightcurveStruct.dataname = None
        CustomLightcurveStruct.plotname = None

    CustomLightcurveStruct.savedata = types.MethodType(savedata_disabled, CustomLightcurveStruct)

    return CustomLightcurveStruct


def CustomImageStruct(source: int = None, pos: list[float] = None) -> CustomImageStruct:
    """CustomImageStruct(source/pos)
    Generates an empty :class:`ImageStruct <AstroToolkit.Data.imagequery.ImageStruct>`, as returned by :ref:`image queries <image-query>`.

    Optionally, the pos/source to which the data structure will refer may be passed. Some attributes will then be automatically filled. If a source or pos is not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: :class:`ImageStruct <AstroToolkit.Data.imagequery.ImageStruct>`
    :rtype: class

    |

    """

    from .Data.imagequery import ImageStruct as BaseImageStruct

    corrected_inputs = check_inputs({"source": [source, int], "pos": [pos, list]}, "customimagestruct")
    source, pos = corrected_inputs

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia").data
        pos = [gaia_data["ra2000"][0], gaia_data["dec2000"][0]]
        identifier = identifier_from_pos(pos)
    elif pos:
        identifier = identifier_from_pos(pos)
    else:
        identifier = None

    CustomImageStruct = BaseImageStruct(survey=None, source=source, pos=pos, identifier=identifier, data=None)

    if source or pos:
        fname = name_file(CustomImageStruct)
        CustomImageStruct.dataname = fname
        generate_plotname(CustomImageStruct)
    else:
        CustomImageStruct.dataname = None
        CustomImageStruct.plotname = None

    CustomImageStruct.savedata = types.MethodType(savedata_disabled, CustomImageStruct)

    return CustomImageStruct


def CustomHrdStruct(sources: list[int] = None) -> CustomHrdStruct:
    """CustomHrdStruct(source/pos)
    Generates an empty :class:`HrdStruct <AstroToolkit.Data.hrdquery.HrdStruct>`, as returned by :ref:`hrd queries <hrd-query>`.

    Optionally, the list of sources to which the data structure will refer may be passed. Some attributes will then be automatically filled. If these are not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param sources: Gaia source IDs
    :type sources: list<int>

    :return: :class:`HrdStruct <AstroToolkit.Data.hrdquery.HrdStruct>`
    :rtype: class

    |

    """
    from .Data.hrdquery import HrdStruct as BaseHrdStruct

    corrected_inputs = check_inputs({"sources": [sources, list]}, "customhrdstruct")
    sources = corrected_inputs[0]

    CustomHrdStruct = BaseHrdStruct(sources=sources, data=None)

    if sources:
        identifiers = []
        for source in sources:
            gaia_data = query(kind="data", source=source, survey="gaia")
            identifiers.append(gaia_data.identifier)
        CustomHrdStruct.identifiers = identifiers
        CustomHrdStruct.dataname = name_file(CustomHrdStruct)
        generate_plotname(CustomHrdStruct)
    else:
        CustomHrdStruct.dataname = None
        CustomHrdStruct.plotname = None
        CustomHrdStruct.identifiers = None

    CustomHrdStruct.savedata = types.MethodType(savedata_disabled, CustomHrdStruct)

    return CustomHrdStruct


def CustomSedStruct(source: int = None, pos: list[float] = None) -> CustomSedStruct:
    """CustomSedStruct(source/pos)
    Generates an empty :class:`SedStruct <AstroToolkit.Data.sedquery.SedStruct>`, as returned by :ref:`sed queries <sed-query>`.

    Optionally, the pos/source to which the data structure will refer may be passed. Some attributes will then be automatically filled. If a source or pos is not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: :class:`SedStruct <AstroToolkit.Data.sedstruct.SedStruct>`
    :rtype: class

    |

    """
    from .Data.sedquery import SedStruct as BaseSedStruct

    corrected_inputs = check_inputs({"source": [source, int], "pos": [pos, list]}, "customsedstruct")
    source, pos = corrected_inputs

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia").data
        pos = [gaia_data["ra2000"][0], gaia_data["dec2000"][0]]
        identifier = identifier_from_pos(pos)
    elif pos:
        identifier = identifier_from_pos(pos)
    else:
        identifier = None

    CustomSedStruct = BaseSedStruct(source=source, pos=pos, identifier=identifier, data=None)

    if source or pos:
        fname = name_file(CustomSedStruct)
        CustomSedStruct.dataname = fname
        generate_plotname(CustomSedStruct)
    else:
        CustomSedStruct.dataname = None
        CustomSedStruct.plotname = None

    CustomSedStruct.savedata = types.MethodType(savedata_disabled, CustomSedStruct)

    return CustomSedStruct


def CustomSpectrumStruct(source: int = None, pos: list[float] = None):
    """CustomSpectrumStruct(source/pos)
    Generates an empty :class:`SpectrumStruct <AstroToolkit.Data.spectrumquery.SpectrumStruct>`, as returned by `spectrum queries <spectrum-query>`.

    Optionally, the pos/source to which the data structure will refer may be passed. Some attributes will then be automatically filled. If a source or pos is not provided, the returned structure will be entirely empty (i.e. all attributes will be set to None).

    :param int source: Target GAIA DR3 Source ID
    :param list<float> pos: Position [right ascension, declination] in degrees

    :return: :class:`SpectrumStruct <AstroToolkit.Data.spectrumquery.SpectrumStruct>`
    :rtype: class

    |

    """
    from .Data.spectrumquery import SpectrumStruct as BaseSpectrumStruct

    corrected_inputs = check_inputs({"source": [source, int], "pos": [pos, list]}, "customspectrumstruct")
    source, pos = corrected_inputs

    if source:
        gaia_data = query(kind="data", source=source, survey="gaia").data
        pos = [gaia_data["ra2000"][0], gaia_data["dec2000"][0]]
        identifier = identifier_from_pos(pos)
    elif pos:
        identifier = identifier_from_pos(pos)
    else:
        identifier = None

    CustomSpectrumStruct = BaseSpectrumStruct(survey=None, source=source, pos=pos, identifier=identifier, data=None)

    CustomSpectrumStruct.savedata = types.MethodType(savedata_disabled, CustomSpectrumStruct)

    from .FileHandling.file_naming import generate_plotname

    if source or pos:
        fname = name_file(CustomSpectrumStruct)
        CustomSpectrumStruct.dataname = fname
        generate_plotname(CustomSpectrumStruct)
    else:
        CustomSpectrumStruct.dataname = None
        CustomSpectrumStruct.plotname = None

    return CustomSpectrumStruct
