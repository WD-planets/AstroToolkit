"""
#################
Manipulating Data
#################
In previous tutorials, we changed the way that we plotted the returned data containers - but we can also manipulate the containers themselves.

Container Methods
=================
The set of available methods depends on the kind of data, but all methods can be utilised via the :meth:`~ATK.Models.DataSet.apply` method of a :class:`~ATK.Models.DataSet`. To showcase this, we will start by fetching some TESS data for SU UMa, the prototype star of a subclass of cataclysmic variables:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=1091051096255456384, survey="tess", path="example_lightcurve_2.fits", split=True)
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
# We can see that TESS has sporadically observed SU UMa over the last ~4 years. But what if we wanted to focus on one of these observation periods? We can crop the light curve by applying the :meth:`~ATK.Models.Lightcurve.crop` method, which truncates all array-like attributes of a class:`~ATK.Models.Lightcurve` to match a given MJD range:

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
# This modifies the :class:`Lightcurves <ATK.Models.Lightcurve>` in-place, but we can also operate on a copy of the :class:`~ATK.Models.DataSet`, leaving the original unmodified. We do this by passing ``inplace=False`` to :meth:`~ATK.Models.DataSet.apply`:

cropped_data = asassn_query.apply("crop", min=60310, max=60340, inplace=False)

# %%
# |
#
# Binning a Light Curve
# ---------------------
# We can also bin a light curve, which combines all array-like attributes of a :class:`~ATK.Models.Lightcurve` into a requested number of ``bins``:

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
# |
#
# or into bins of a given ``size``:

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
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
