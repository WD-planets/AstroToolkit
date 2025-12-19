import astropy.units as u
from astropy.coordinates import SkyCoord

from ATK.Tools import Target, query, read

SOURCE = 2552928187080872832
POSITION = SkyCoord(ra=12.297, dec=5.376, unit=u.deg, frame="icrs")

target = Target.from_id(SOURCE, "gaia")
# target = Target.from_pos(POSITION)

data = query(kind="vizier", target=target, survey="galex", radius=5)

path = data.save()
data = read(path)

data.show(show_all_types=True)
