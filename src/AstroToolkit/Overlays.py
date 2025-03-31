from .Configuration.overlays import OverlayStruct

overlays = OverlayStruct()


def addOverlay(survey, ra_name, dec_name, id_name, mag_names=None):
    """addOverlay(survey,ra_name,dec_name,id_name,mag_names)
    Adds an overlay definition to ATK for use in image queries and plotting.

    :param survey: name of a catalogue alias, e.g. 'allwise'
    :type name: str
    :param ra_name: name of right ascension column in Vizier
    :type ra_name: str
    :param dec_name: name of declination column in Vizier
    :type dec_name: str
    :param id_name: name of survey ID column in Vizier
    :type id_name: str
    :param mag_names: names of magnitude columns in Vizier
    :type mag_names: list<str>, optional

    :return: None

    |

    """
    overlays.add_overlay(survey, ra_name, dec_name, id_name, mag_names)
    return None


def delOverlay(kind, survey):
    """delOverlay(kind,survey)
    Deletes an existing overlay definition.

    :param kind: overlay kind, from: scaled_detection, detection
    :type kind: str
    :param survey: overlay survey
    :type survey: str

    :return: None

    |

    """
    overlays.del_overlay(kind, survey)

    print(f"Deleted existing {kind} overlay definition for survey `{survey}`.")
    return None


def resetOverlays():
    """resetOverlays()
    Resets the list of ATK overlay definitions to its default state.

    :return: None

    |

    """
    overlays.reset_overlays()
    print("Resetting ATKOverlays.yaml to default values...\n")
    return None


def openOverlays():
    """openOverlays()
    Opens the list of ATK overlay definitions in the default text editor.

    :return: None

    |

    """
    path = overlays.overlay_file

    import platform
    import subprocess

    if platform.system().lower() in ["posix", "linux"]:
        subprocess.run(["chmod", "+x", str(path)])
        subprocess.run(["xdg-open", str(path)])
    else:
        import webbrowser

        webbrowser.open(path)

    return None
