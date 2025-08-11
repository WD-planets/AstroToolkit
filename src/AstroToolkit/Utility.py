def openFileDialogue():
    import os

    from PyQt5.QtWidgets import QApplication, QFileDialog

    app = QApplication([])
    fname, _ = QFileDialog.getOpenFileName(None, "Open File", str(os.getcwd()), "ATK Data File (*.fits)")

    return fname


def HJDtoMJD(hjd, pos):
    import astropy.units as u
    from astropy.coordinates import SkyCoord, get_sun
    from astropy.time import Time

    def get_mjd(hjd, pos):
        t_hjd = Time(hjd, format="jd", scale="utc")
        sun_position = get_sun(t_hjd)
        target = SkyCoord(ra=pos[0] * u.deg, dec=pos[1] * u.deg, frame="icrs")
        heliocentric_correction = sun_position.cartesian.dot(target.cartesian).to(u.au).value / 173.144632674240
        jd = t_hjd.jd - heliocentric_correction

        return Time(jd, format="jd", scale="utc").mjd

    if isinstance(hjd, list):
        calculated_mjds = []
        for entry in hjd:
            calculated_mjds.append(get_mjd(entry, pos))
        return calculated_mjds
    elif isinstance(hjd, (float, int)):
        return get_mjd(hjd, pos)


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
