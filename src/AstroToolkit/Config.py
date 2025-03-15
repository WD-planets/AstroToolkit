"""
Most defaults in ATK can be changed via a config file. This module allows for the configuration file to be viewed and edited.
"""

from .Configuration.baseconfig import ConfigStruct
from .Input.input_validation import check_inputs

config = ConfigStruct()
config.read_config()


def editConfig(key: str, value: str) -> None:
    """editConfig(key, value)
    Edits config values. See :ref:`Config Keys` for a list and description of available keys.

    :param key: config key
    :type key: str
    :param value: value to assign to this config key
    :type value: str

    :return: None

    |

    """

    corrected_inputs = check_inputs(
        {"key": [key, str], "value": [value, str]}, "editconfig"
    )
    key, value = corrected_inputs

    print("Written change to ATKConfig.ini. New Values:\n")
    config.edit_config(key, value)

    return None


def openConfig() -> None:
    """openConfig()
    Opens the config in the default text editor. See :ref:`Config Keys` for a list and description of available keys.

    :return: None

    |

    """
    config = ConfigStruct()
    path = config.config_file

    import platform
    import subprocess

    if platform.system().lower() in ["posix", "linux"]:
        subprocess.run(["chmod", "+x", str(path)])
        subprocess.run(["xdg-open", str(path)])
    else:
        import webbrowser

        webbrowser.open(path)

    return None


def showConfig() -> None:
    """showConfig()
    Prints the current config file to stdout.

    :return: None

    |

    """
    print("Current ATKConfig.ini values:\n")
    config.read_config()
    config.output_config()

    return None


def resetConfig() -> None:
    """resetConfig()
    Resets the config to default values. A list of available keys and their default values can be found in :ref:`Config Keys`.

    :return: None

    |

    """
    print("Resetting ATKConfig.ini to default values...\n")
    config.set_default_config()
    config.write_config()

    return None
