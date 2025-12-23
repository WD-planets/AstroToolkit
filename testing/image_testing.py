import astropy.units as u
from astropy.coordinates import SkyCoord

import ATK.Tools as ATK

SOURCE = 587316166180416640
SOURCE = 2552928187080872832
POS = SkyCoord(ra=141.18533 * u.deg, dec=8.03083 * u.deg, frame="icrs")
GAL_POS = SkyCoord(l=224.26315 * u.deg, b=37.60412 * u.deg, frame="galactic")

target = ATK.Target.from_id(SOURCE)

data = ATK.query("image", survey="panstarrs", target=GAL_POS, size=120, band="g", overlays=["gaia", "galex"], disable_corrections=True)

# data = ATK.query("image", survey="dss1", target=SOURCE, size=300, band="blue", overlays=["galex"])

if data.data:
    img = data.data[0]
    img.show()

    print("\n\n---------------------------------\n\n")

    data.show()
    data.plot(relative_axes=False)
    data.open()
