"""
#########################
Light Curve Customisation
#########################
When plotting light curves, :meth:`~ATK.Models.DataSet.plot` - and hence :meth:`~ATK.Models.DataSet.open` - accepts the following additional arguments.

|

Specifying Bands
================
If only certain photometric bands are of interest, the ``bands`` parameter can be passed as a list of bands that should be included in plotting:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits.gz")
# sphinx_gallery_start_ignore
asassn_query.plot(bands=["g"])
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open(bands=["g"])
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    ATK supports the following lightcurve surveys and bands:
#
#    .. include:: supported_lightcurve_surveys.rst
#
# |
# 
# Setting Band Colours
# ====================
# The colours of each band can also be chosen:

# sphinx_gallery_start_ignore
asassn_query.plot(bands=["v","g"], colours=["orange","blue"])
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open(bands=["v","g"], colours=["orange","blue"])
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. note::
#    The following colours are supported: ``"green"``, ``"red"``, ``"blue"``, ``"orange"``, ``"purple"``, ``"black"``

# %%
# |
#
# Choosing a Colour Map
# =====================
# By default, the colour map scales with distance from the mean brightness of the photometry in each band. This can be disabled by setting ``cmap = "flat"`` (default = ``"mean"``):

# sphinx_gallery_start_ignore
asassn_query.plot(cmap="flat")
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open(cmap="flat")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Using Non-Reduced MJD
# =====================
# By default, the x-axis is reduced to show the time since the earliest observation in each set of light curves (i.e. the minimum MJD is subtracted from all bands). To instead plot against MJD, pass ``time_format = "original"`` (default = ``"reduced"``):

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
# |
# |
#
# .. rubric:: Download this Tutorial
