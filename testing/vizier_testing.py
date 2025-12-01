import astropy.units as u
from astropy.coordinates import SkyCoord

from ATK.structures.data_containers.definitions import Image
from ATK.Tools import query

SOURCE = 2552928187080872832
POSITION = SkyCoord(ra=12.291, dec=05.389, unit=u.deg, frame="icrs")
GAL_POSITION = SkyCoord(l=121.880, b=-57.478, unit=u.deg, frame="galactic")

data = query(kind="vizier", target=POSITION, survey="gaia", radius=100)

additional_data = True
if additional_data:
    img = Image()
    img2 = Image()

    data.data = [img, img2]


data.show()

# data = query(kind="vizier", target=SOURCE, survey="galex")
# data.show()

"""
data = query(kind="vizier", target=POSITION, survey="gaia")
data.show()

data = query(kind="vizier", target=POSITION, survey="galex")
data.show()

data = query(kind="vizier", target=GAL_POSITION, survey="gaia")
data.show()

data = query(kind="vizier", target=GAL_POSITION, survey="galex")
data.show()
"""
