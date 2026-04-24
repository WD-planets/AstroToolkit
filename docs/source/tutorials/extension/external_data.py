"""
#######################
Importing External Data
#######################
Generating a Dataset
====================
ATK also supports the use of external data. The first step is to create an empty :class:`~ATK.Models.DataSet`. This is done with :meth:`~ATK.Models.DataSet.from_target`, which requires any valid targeting for initialisation (see :doc:`here <../getting_started/data_query>`). For example, using a Gaia Source ID:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
import pandas as pd
import numpy as np
from astropy.io import fits
# sphinx_gallery_end_ignore

from ATK.Models import DataSet

target = 4223502720986764672

dataset = DataSet.from_target(target)

# %%
# |
#
# The next step is to generate a data container to add to the :class:`~ATK.Models.DataSet`. For this tutorial, an externally-acquired spectrum will be imported into ATK. The easiest way to do this is to generate a :class:`~ATK.Models.Spectrum` using one of its IO methods: :meth:`~ATK.Models.Spectrum.from_table` and :meth:`~ATK.Models.Spectrum.from_dataframe`. In this tutorial, a :class:`~pandas.DataFrame` will be used.
# 
# The :class:`~pandas.DataFrame` (or :class:`~astropy.table.Table` if using :meth:`~ATK.Models.Spectrum.from_table`) must contain a column with the same name as each array-like attribute of the required data container. For a :class:`~ATK.Models.Spectrum`, the :class:`~pandas.DataFrame` must have a ``wavelength`` columns and a ``flux`` column. 
# 
# .. note:: 
# 
#    When generating data containers with :meth:`~ATK.base.DataFrameIOMixin.from_dataframe`, all columns are assumed to be in specific units - details for which can be found in the container's documentation (e.g. :class:`here <ATK.Models.Spectrum>` for spectra). While units could be manually changed after initialisation, it may instead be preferable to work with a :class:`~astropy.table.Table`, as :meth:`~ATK.Models.Spectrum.from_table` preserves all units of the input :class:`~astropy.table.Table`.
# 
# |
# 
# Since the format of external data can vary greatly, the process of generating a :class:`~pandas.DataFrame` will not be covered here, but this is what the :class:`~pandas.DataFrame` should look like:
# 
# .. code-block:: python
# 
#    df = ...

# sphinx_gallery_start_ignore
hdul = fits.open("external_spec.fits")
header = hdul[0].header
flux = hdul[0].data

crval1 = header["CRVAL1"]
crpix1 = header["CRPIX1"]
cdelt1 = header["CDELT1"]
naxis1 = header["NAXIS1"]
mjd = header["MJD-OBS"]

pixel_indices = np.arange(naxis1)
wavelength = crval1 + (pixel_indices + 1 - crpix1) * cdelt1

df = pd.DataFrame({"wavelength": wavelength * 10, "flux": flux * 10**17})
# sphinx_gallery_end_ignore

print(df)

# %%
# |
# 
# Any other required attributes are passed as keyword arguments to :meth:`~ATK.Models.Spectrum.from_dataframe` (or :meth:`~ATK.Models.Spectrum.from_table`). A :class:`~ATK.Models.Spectrum` only requires that we specify a ``survey``. 
# 
# .. note::
# 
#    For information on the keyword arguments that are required by a data container, see its documentation for the IO method that is being used (for this tutorial, that would be :meth:`~ATK.Models.Spectrum.from_dataframe`).
# 
# |
#
# The following will generate a :class:`~ATK.Models.Spectrum` which contains the external data:

from ATK.Models import Spectrum

spec = Spectrum.from_dataframe(target, data=df, survey="XSHOOTER")
spec.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# |
# 
# The final step is to add the :class:`~ATK.Models.Spectrum` to the :class:`~ATK.Models.DataSet`, which can be done with :meth:`~ATK.Models.DataSet.add`:

dataset.add(spec)
dataset.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# |
#
# The newly constructed :class:`~ATK.Models.DataSet` can then be used throughout ATK as if it were retrieved internally:

dataset.apply("crop", min=5700, inplace=True)
# sphinx_gallery_start_ignore
dataset.plot(fit=True, prom=2, smooth=5)
figure = format_plot(dataset.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
dataset.open(fit=True, prom=2, smooth=5)
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Another Example
# ===============
# Below is a full example of the process for importing and utilising an externally-acquired light curve (excluding loading of the data):
# 
# .. code-block:: python
# 
#    df = ...

# sphinx_gallery_start_ignore
def hjd_to_mjd(hjd):
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import get_sun
    from astropy.time import Time

    hjd = np.asarray(hjd)

    # Create Time object (vectorized)
    hjd_time = Time(hjd, format="jd")

    # Get Sun position (vectorized)
    sun_pos = get_sun(hjd_time)

    # Light travel time correction (in days)
    heliocentric_correction = (sun_pos.distance.to(u.au).value / 1731.456) * u.day

    # Convert HJD → JD
    jd_time = hjd_time - heliocentric_correction

    # Return MJD directly
    return jd_time.mjd

df = pd.read_csv(
    "external_lightcurve.dat", names=["hjd", "zero", "mag", "mag_err"], delimiter=" "
)
df["mjd"] = hjd_to_mjd(df["hjd"])
# sphinx_gallery_end_ignore
from ATK.Models import Lightcurve

target = 6050296829033196032 

dataset = DataSet.from_target(target)
lc = Lightcurve.from_dataframe(target, data=df, survey="TNT", band="KG5")
dataset.add(lc)

folded = dataset.apply("fold", fmin=0, fmax=10, samples=10000, inplace=False)
folded.apply("bin", bins=200)
# sphinx_gallery_start_ignore
folded.plot()
figure = format_plot(folded.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
