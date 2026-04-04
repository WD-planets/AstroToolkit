"""
#################
Manipulating Data
#################
:class:`Lightcurves <ATK.Models.Lightcurve>` offer the first examples of data methods - these are :class:`~ATK.Models.DataSet` methods that manipulate the containers themselves.

Data Methods
============
The set of available methods depends on the kind of data that is being stored in the :class:`~ATK.Models.DataSet`, but all methods can be utilised via the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.apply` method.

To showcase this, some data will first be fetched from TESS for SU UMa, the prototype star of a subclass of cataclysmic variables:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=1091051096255456384, survey="tess", path="example_lightcurve_2.fits.gz")
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
# According to the above, TESS has sporadically observed SU UMa over the last 4-5 years. To focus on one of these periods, the light curve can be cropped by applying :meth:`~ATK.Models.Lightcurve.crop`, which truncates all array-like attributes of the :class:`~ATK.Models.Lightcurve` to match a chosen MJD range:

asassn_query.apply("crop", min=60310, max=60340)
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
# 
# |
# 
# **By default,** :meth:`~ATK.Models.DataSet.apply` **modifies the stored containers in-place**. **To instead operate on a copy of the** :class:`~ATK.Models.DataSet` **- thereby leaving the original unmodified - pass** ``inplace=False`` **to** :meth:`~ATK.Models.DataSet.apply`:

cropped_data = asassn_query.apply("crop", min=60310, max=60340, inplace=False)

# %%
# |
#
# Binning a Light Curve
# ---------------------
# :class:`Lightcurves <ATK.Models.Lightcurve>` can also be binned, combining all array-like attributes into a requested number of ``bins`` or a given bin ``size``:

binned_data = asassn_query.apply("bin", bins=1000, inplace=False)
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

binned_data = asassn_query.apply("bin", size=1*u.hour, inplace=False)
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
# :class:`Lightcurves <ATK.Models.Lightcurve>` can be sigma clipped by applying :meth:`~ATK.Models.Lightcurve.clip`, which sigma clips all array-like attributes to eliminate data points where the brightness is outside a given sigma range:

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
