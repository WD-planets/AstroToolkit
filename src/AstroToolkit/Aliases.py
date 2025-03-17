"""
ATK supports the creation of aliases to any Vizier catalogue, which are set via an aliases file. This module allows the aliases file to be viewed and edited.
"""

from .Configuration.catalogue_setup import CatalogueStruct
from .Input.input_validation import check_inputs

catalogues = CatalogueStruct()
catalogues.get_catalogues()


def addAlias(name: str, id: str) -> None:
    """addAlias(name,id)
    Adds a Vizier catalogue alias to ATK for use in data queries.

    :param name: alias name, e.g. 'allwise'
    :type name: str
    :param id: Vizier catalogue ID, e.g. 'II/328/allwise'
    :type id: str

    :return: None

    |

    """

    corrected_inputs = check_inputs({"name": [name, str], "id": [id, str]}, "addAlias")
    name, id = corrected_inputs

    catalogues.add_catalogue(name, id)
    print(f"Added alias for {id} with label {name}.")

    return None


def delAlias(name: str) -> None:
    """delAlias(name)
    Deletes an existing catalogue alias.

    :param name: alias name to be deleted
    :type name: str

    :return: None

    |

    """

    corrected_inputs = check_inputs({"name": [name, str]}, "delAlias")
    name = corrected_inputs[0]

    catalogues.delete_catalogue(name)
    print(f"Deleted alias with label {name}.")


def openAliases() -> None:
    """openAliases()
    Opens the catalogue alias list in the default text editor.

    :return: None

    |

    """
    catalogues = CatalogueStruct()
    path = catalogues.catalogue_file

    import platform
    import subprocess

    if platform.system().lower() in ["posix", "linux"]:
        subprocess.run(["chmod", "+x", str(path)])
        subprocess.run(["xdg-open", str(path)])
    else:
        import webbrowser

        webbrowser.open(path)

    return None


def resetAliases() -> None:
    """resetAliases()
    Resets the catalogue alias list (i.e. only keeps default ATK data surveys).

    :return: None

    |

    """
    print("Resetting ATKAliases.ini to default values...\n")
    catalogues.default_setup()

    return None


def showAliases() -> None:
    """showAliases()
    Prints the current catalogue alias list to stdout.

    :return: None

    |

    """
    print("Current ATKAliases.ini values:\n")
    catalogues.get_catalogues()
    catalogues.output_catalogues()

    return None
