"""
Any epochs used in the correction of proper motion throughout ATK can be set via an epochs file. This module allows the epochs file to be viewed and edited.
"""

from .Configuration.epochs import EpochStruct
from .Input.input_validation import check_inputs

epoch_struct = EpochStruct()


def setEpoch(section: str, survey: str, epoch: list) -> None:
    """setEpoch(section,survey,epoch)
    Sets the epoch of a supported data, lightcurve or spectrum survey, or a Vizier catalogue alias.

    :param section: section name from: default_data_surveys,
                    additional_data_surveys, lightcurve_surveys,
                    spectrum_surveys
    :type section: str
    :param survey: survey for which an epoch definition should be set
    :type survey: str
    :param epoch: survey epoch in format: [year,month]
    :type epoch: list<int>

    :return: None

    |

    """

    corrected_inputs = check_inputs(
        {"section": [section, str], "survey": [survey, str], "id": [epoch, list]}, "setEpoch"
    )
    section, survey, epoch = corrected_inputs

    epoch_struct.set_epoch(section=section, survey=survey, epoch=epoch)

    return None


def delEpoch(alias: str) -> None:
    """delEpoch(name)
    Deletes an existing epoch definition for a Vizier catalogue alias.

    :param alias: alias name for which an epoch definition should be deleted
    :type alias: str

    :return: None

    |

    """

    corrected_inputs = check_inputs({"survey": [alias, str]}, "delEpoch")
    alias = corrected_inputs[0]

    print(f"Deleted existing epoch definition for alias '{alias}'.")

    epoch_struct.delete_epoch(alias)


def openEpochs() -> None:
    """openEpochs()
    Opens the epoch definition file in the default text editor.

    :return: None

    |

    """
    catalogues = epoch_struct
    path = catalogues.epoch_file

    import platform
    import subprocess

    if platform.system().lower() in ["posix", "linux"]:
        subprocess.run(["chmod", "+x", str(path)])
        subprocess.run(["xdg-open", str(path)])
    else:
        import webbrowser

        webbrowser.open(path)

    return None


def resetEpochs() -> None:
    """resetEpochs()
    Resets the catalogue alias list.

    :return: None

    |

    """
    print("Resetting ATKEpochs.ini to default values...\n")
    epoch_struct.default_setup()

    return None


def showEpochs() -> None:
    """showEpochs()
    Prints the current epoch defintions to stdout.

    :return: None

    |

    """
    print("Current ATKEpochs.ini values:\n")
    epoch_struct.output_epochs()

    return None
