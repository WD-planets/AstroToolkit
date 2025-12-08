import numpy as np

from ATK.structures.definitions import Lightcurve

mjd = np.array([1, 2, 3])
flux = np.array([1, 2, 3])
flux_err = np.array([1, 2, 3])

lc = Lightcurve("ztf", "g", mjd=mjd, flux=flux, flux_err=flux_err)

print(lc.__dict__)
