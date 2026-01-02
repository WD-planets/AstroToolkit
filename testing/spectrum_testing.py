import astropy.units as u
from astropy.coordinates import SkyCoord

import ATK.Tools as ATK

POS = SkyCoord(ra=250.423475, dec=36.461319, unit=(u.deg, u.deg))

SOURCE = 587316166180416640

desi_data = ATK.query("spectrum", survey="desi", target=SOURCE, radius=100)
desi_data.save("test_desi_spectrum.fits")
rec_data = ATK.read("test_desi_spectrum.fits")
rec_data.show()
rec_data.open()

# data = ATK.query("spectrum", survey="sdss", target=SOURCE, radius=100)
# data.data += desi_data.data

# data.show()
# data.open()
