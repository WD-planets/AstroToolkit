import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

from ...structures.definitions import Spectrum, Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS


def query(target: Target, **kwargs: dict):
    radius = kwargs["radius"] * u.arcsec

    # get any SDSS spectra in radius
    try:
        data = SDSS.get_spectra(coordinates=target.coords, radius=radius, timeout=180)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    spectra = []
    for entry in data:
        primary_hdr = entry[0].header
        exposure = primary_hdr["EXPTIME"]
        spec_pos = SkyCoord(ra=primary_hdr["PLUG_RA"] * u.deg, dec=primary_hdr["PLUG_DEC"] * u.deg, frame=primary_hdr["RADECSYS"].lower())

        data = entry[1].data
        wavelength = 10 ** data["loglam"]
        # fluxes are already in e-17 erg/cm^2/s/Ang
        flux = data["flux"]

        spec = Spectrum(
            "sdss",
            wavelength=wavelength,
            flux=flux,
            exposure=exposure,
            position=spec_pos,
            separation=spec_pos.separation(target.coords).to(u.arcsec).value,
        )
        spectra.append(spec)

    return spectra
