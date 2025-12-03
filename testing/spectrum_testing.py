import numpy as np

from ATK.structures.data_containers.definitions import Spectrum

spec = Spectrum()

spec.survey = "sdss"
spec.flux = np.array([1, 2, 3])
spec.wavelength = np.array([1, 2, 3])

df = spec.to_dataframe()
hdu = spec.to_hdu()

print(df)
print(hdu.data)
