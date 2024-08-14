"""
This module allows for the creation of empty ATK data structures, which can be filled with custom data for use throughout ATK. These structures will initialise with all attributes set to None.
"""


def CustomDataStruct():
    """
    Returns an empty :ref:`DataStruct`, as returned by data queries.

    :return: empty ATK DataStruct
    :rtype: class

    |

    """

    from .Data.dataquery import DataStruct as BaseDataStruct

    CustomDataStruct = BaseDataStruct(
        survey=None,
        catalogue=None,
        source=None,
        pos=None,
        identifier=None,
        data=None,
    )

    CustomDataStruct.subkind = None

    return CustomDataStruct


def CustomLightcurveStruct():
    """
    returns an empty :ref:`LightcurveStruct`, as returned by lightcurve queries.

    :return: empty ATK LightcurveStruct
    :rtype: class

    |

    """

    from .Data.lightcurvequery import LightcurveStruct as BaseLightcurveStruct

    CustomLightcurveStruct = BaseLightcurveStruct(
        survey=None, source=None, pos=None, identifier=None, data=None
    )

    return CustomLightcurveStruct


def CustomImageStruct():
    """
    Returns an empty :ref:`ImageStruct`, as returned by image queries.

    :return: empty ATK ImageStruct
    :rtype: class

    |

    """

    from .Data.imagequery import ImageStruct as BaseImageStruct

    CustomImageStruct = BaseImageStruct(survey=None, source=None, pos=None, data=None)

    return CustomImageStruct


def CustomHrdStruct():
    """
    Returns an empty :ref:`HrdStruct`, as returned by hrd queries.

    :return: empty ATK HrdStruct
    :rtype: class

    |

    """
    from .Data.hrdquery import HrdStruct as BaseHrdStruct

    CustomHrdStruct = BaseHrdStruct(sources=None, data=None)

    return CustomHrdStruct


def CustomSedStruct():
    """
    Returns an empty :ref:`SedStruct`, as returned by sed queries.

    :return: empty ATK SedStruct
    :rtype: class

    |

    """
    from .Data.sedquery import SedStruct as BaseSedStruct

    CustomSedStruct = BaseSedStruct(source=None, pos=None, data=None)

    return CustomSedStruct


def CustomSpectrumStruct():
    """
    Returns an empty :ref:`SpectrumStruct`, as returned by spectrum queries.

    :return: empty ATK SpectrumStruct
    :rtype: class

    |

    """
    from .Data.spectrumquery import SpectrumStruct as BaseSpectrumStruct

    CustomSpectrumStruct = BaseSpectrumStruct(
        survey=None, source=None, pos=None, data=None
    )

    return CustomSpectrumStruct
