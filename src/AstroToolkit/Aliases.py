from .Configuration.catalogue_setup import CatalogueStruct

catalogues = CatalogueStruct()
catalogues.get_catalogues()


def addAlias(name: str, id: str) -> None:
    """addAlias(name,id)
    Adds a Vizier catalogue alias to ATK for use in data queries.

    :param name: catalogue name, e.g. 'allwise'
    :type name: str
    :param id: Vizier catalogue ID, e.g. 'II/328/allwise'
    :type id: str

    :return: None

    |

    """

    catalogues.add_catalogue(name, id)
    print(f"Added alias for {id} with label {name}.")

    return None


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
