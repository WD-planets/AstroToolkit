def openFileDialogue():
    import os

    from PyQt5.QtWidgets import QApplication, QFileDialog

    app = QApplication([])
    fname, _ = QFileDialog.getOpenFileName(None, "Open File", str(os.getcwd()), "ATK Data File (*.fits)")

    return fname


def HJDtoMJD(hjd):
    import astropy.units as u
    from astropy.coordinates import get_sun
    from astropy.time import Time

    def get_mjd(hjd):
        hjd_time = Time(entry, format="jd")
        sun_position = get_sun(hjd_time)
        heliocentric_correction = sun_position.distance.to(u.au).value / 1731.456
        jd = hjd_time - heliocentric_correction * u.day

        return (jd - 2400000.5).value

    if isinstance(hjd, list):
        calculated_mjds = []
        for entry in hjd:
            calculated_mjds.append(get_mjd(entry))
        return calculated_mjds
    elif isinstance(hjd, float):
        return get_mjd(hjd)


def getBrightnessType(data):
    brightness_types = []
    for band in data:
        if "mag" in band:
            if "flux_err" in band:
                raise ValueError("Invalid combination of 'mag' and 'flux_err'.")
            brightness_type = "mag"
            brightness_types.append("mag")
        elif "flux" in band:
            if "mag_err" in band:
                raise ValueError("Invalid combination of 'flux' and 'mag_err'.")
            brightness_type = "flux"
            brightness_types.append("flux")
        else:
            raise ValueError("Invalid brightness type, expected 'mag' and 'mag_err' or 'flux' and 'flux_err'.")

    if len(list(dict.fromkeys(brightness_types))) > 1:
        raise ValueError(
            "Inconsistent brightness types among bands. Expected consistent use of 'mag' and 'mag_err' or 'flux' and 'flux_err'."
        )

    return brightness_type
