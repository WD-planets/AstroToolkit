"""
#################
Manipulating Data
#################
:class:`~ATK.Models.Lightcurve` objects offer the first examples of **data methods** - these are methods that manipulate the containers themselves.

Data Methods
============
The set of available methods depends on the ``kind`` of data that is being stored in the :class:`~ATK.Models.DataSet`, but all methods can be utilised via the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.apply` method.

|

To showcase this, some data will first be fetched from ZTF for Al Com - a WZ Sge cataclysmic variable:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=3932951035266314496 , survey="ztf", path="example_lightcurve_2.fits.gz")
asassn_query.show(show_all=True)
# sphinx_gallery_start_ignore
asassn_query.plot(time_format="original")
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open(time_format="original")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Cropping a Light Curve
# ----------------------
# The returned :class:`~ATK.Models.Lightcurve` can be cropped by applying its :meth:`~ATK.Models.Lightcurve.crop` **data method**. This truncates all array-like attributes of the :class:`~ATK.Models.Lightcurve` to match a chosen MJD range, in this case focusing on the large outburst:

cropped_data = asassn_query.apply("crop", cmin=58550, cmax=58660, inplace=False)
# sphinx_gallery_start_ignore
cropped_data.plot(time_format="original")
figure = format_plot(cropped_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
cropped_data.open(time_format="original")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# 
# |
# 
# .. note:: 
# 
#    **By default,** :meth:`~ATK.Models.DataSet.apply` **modifies the stored data containers in-place**. To instead operate on a copy of the :class:`~ATK.Models.DataSet` - thereby leaving the original unmodified - pass ``inplace=False`` to :meth:`~ATK.Models.DataSet.apply` (as above).
#
# |
#
# Binning a Light Curve
# ---------------------
# The :class:`~ATK.Models.Lightcurve` can also be binned, combining all array-like attributes into a requested number of ``bins`` or a given bin ``size``:

binned_data = asassn_query.apply("bin", bins=200, inplace=False)
# sphinx_gallery_start_ignore
binned_data.plot(time_format="original")
figure = format_plot(binned_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
binned_data.open(time_format="original")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%

import astropy.units as u

binned_data = asassn_query.apply("bin", size=10*u.day, inplace=False)
# sphinx_gallery_start_ignore
binned_data.plot(time_format="original")
figure = format_plot(binned_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
binned_data.open(time_format="original")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# Sigma-Clipping a Light Curve
# ============================
# The :class:`~ATK.Models.Lightcurve` can be sigma clipped by applying :meth:`~ATK.Models.Lightcurve.clip`, which clips all array-like attributes to eliminate data points where the light curve's brightness is outside a given number of standard deviation from the median:

clipped_data = binned_data.apply("clip", sigma=2, inplace=False)
# sphinx_gallery_start_ignore
clipped_data.plot(time_format="original")
figure = format_plot(clipped_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
clipped_data.open(time_format="original")
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
