"""
Most defaults in ATK can be changed via a config file. This module allows for the configuration file to be viewed and edited.
"""

from .Configuration.baseconfig import ConfigStruct

config = ConfigStruct()
config.read_config()


def editconfig(key, value):
    """
    Allows for config values to be edited. See :ref:`Config Keys` for a description of available keys.

    :param str key: config key
    :param str value: value to assign to this config key

    :return: None

    |

    """
    from .Input.input_validation import check_inputs

    inputs = check_inputs({"key": key, "value": value}, "edit")
    key, value = inputs["key"], inputs["value"]

    print("Written change to ATKConfig.ini. New Values:\n")
    config.edit_config(key, value)


def openconfig():
    """
    Opens the config in the default text editor. See :ref:`Config Keys` for a description of available keys.

    :return: None

    |

    """
    config = ConfigStruct()
    path = config.config_file

    import webbrowser

    webbrowser.open(path)


def outputconfig():
    """
    Prints the current config file to stdout.

    :return: None

    |

    """
    print("Current ATKConfig.ini values:\n")
    config.read_config()
    config.output_config()


def resetconfig():
    """
    Resets the config to default values.

    :return: None

    |

    """
    print("Resetting ATKConfig.ini to default values...\n")
    config.set_default_config()
    config.write_config()
