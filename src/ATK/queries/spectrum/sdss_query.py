import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

from ...structures.definitions import Spectrum, Target
from ...utilities.defaults import CONNECTION_ERRORS, RETURNS
from ...utilities.misc import angle_to_quantity


def query(target: Target, **kwargs: dict):
    """
    Performs an SDSS spectrum query
    """

    radius = kwargs["radius"]

    # get any SDSS spectra in radius
    try:
        data = SDSS.get_spectra(coordinates=target.coords, radius=radius, timeout=180)
    except CONNECTION_ERRORS:
        return RETURNS.EXCEPTION

    if not data:
        return RETURNS.NULL

    spectra = []
    for entry in data:
        # primary header has ancillary information
        primary_hdr = entry[0].header
        exposure = primary_hdr["EXPTIME"]
        spec_pos = SkyCoord(
            ra=primary_hdr["PLUG_RA"] * u.deg,
            dec=primary_hdr["PLUG_DEC"] * u.deg,
            frame=primary_hdr["RADECSYS"].lower(),
        )

        # hdu 1 has spectrum
        data = entry[1].data
        wavelength = 10 ** data["loglam"]

        # fluxes are already in e-17 erg/cm^2/s/Ang
        flux = data["flux"]

        # create Spectrum object for each returned spectrum
        spec = Spectrum(
            survey="sdss",
            wavelength=wavelength * u.Unit("Angstrom"),
            flux=flux * u.Unit("1e-17 erg cm-2 s-1 Angstrom-1"),
            exposure=exposure * u.s,
            search_pos=target.coords,
            # this would normally be an angle
            separation=angle_to_quantity(spec_pos.separation(target.coords), kwargs["radius"].unit),
        )
        spectra.append(spec)

    return spectra
