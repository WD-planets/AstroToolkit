import astropy.units as u
from astropy.coordinates import SkyCoord

import ATK.Tools as ATK

POS = SkyCoord(ra=250.423475, dec=36.461319, unit=(u.deg, u.deg))

SOURCE = 587316166180416640

data = ATK.query("spectrum", survey="sdss", target=SOURCE, radius=100)
data.show()
data.open()
