from .Configuration.overlays import OverlayStruct

overlays = OverlayStruct()


def addOverlay(survey, ra_name, dec_name, id_name, mag_names=None):
    overlays.add_overlay(survey, ra_name, dec_name, id_name, mag_names)
    return None


def delOverlay(kind, survey):
    overlays.del_overlay(kind, survey)

    print(f"Deleted existing {kind} overlay definition for survey `{survey}`.")
    return None


def resetOverlays():
    overlays.reset_overlays()
    print("Resetting ATKOverlays.yaml to default values...\n")
    return None


def openOverlays():
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
