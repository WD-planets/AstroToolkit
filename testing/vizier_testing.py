import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from ATK.structures.definitions import Spectrum
from ATK.Tools import query, read

SOURCE = 2552928187080872832
POSITION = SkyCoord(ra=12.291, dec=05.389, unit=u.deg, frame="icrs")
GAL_POSITION = SkyCoord(l=121.880, b=-57.478, unit=u.deg, frame="galactic")

# data = query(kind="vizier", target=POSITION, survey="gaia", radius=3)
# data.show()

data = query(kind="vizier", target=SOURCE, survey="galex")
spec = Spectrum()
spec.survey = "test"
spec.wavelength = np.array([1, 2, 3])
spec.flux = np.array([1, 2, 3])
# data.data = [spec]
# data.show()
path = data.save()

data = read(path)


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
