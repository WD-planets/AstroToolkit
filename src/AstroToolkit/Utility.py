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
