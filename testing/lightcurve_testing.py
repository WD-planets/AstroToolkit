import astropy.units as u
from astropy.coordinates import SkyCoord

import ATK.Tools as ATK

SOURCE = 587316166180416640
coord = SkyCoord(141.18533 * u.deg, 8.03083 * u.deg, frame="icrs")

ztf_data = ATK.query("lightcurve", survey="ztf", target=SOURCE)
gaia_data = ATK.query("lightcurve", survey="gaia", target=SOURCE)

ztf_data.data += gaia_data.data

ztf_data.plot()
ztf_data.show()
ztf_data.open()
