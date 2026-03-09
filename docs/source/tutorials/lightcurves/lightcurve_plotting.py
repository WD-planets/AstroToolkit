"""
#########################
Light Curve Customisation
#########################

When plotting light curves, :meth:`~ATK.Models.DataSet.plot` accepts a few additional arguments.

|

Selecting Bands
===============
If we only want to plot specific bands of the :class:`Lightcurves <ATK.Models.Lightcurve>` returned by our query, we can pass the ``bands`` parameter to :meth:`~ATK.Models.DataSet.plot` or :meth:`~ATK.Models.DataSet.open`:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits")
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
# |
# 
# Setting Band Colours
# ====================
# We can also choose the colour of each band:

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
#    The following colours are supported:
#    
#    - ``"green"``
#    - ``"red"``
#    - ``"blue"``
#    - ``"orange"``
#    - ``"purple"``
#    - ``"black"``

# %%
# |
#
# Setting the Colour Map
# ======================
# By default the colour map scales with distance from the mean brightness of the photometry in a given band. This can be disabled by setting ``cmap = "flat"`` (default = ``"mean"``):

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
# By default, the x-axis is reduced to show the time since the earliest observation (i.e. the minimum MJD is subtracted). If we instead want to plot the MJD as-is, we can pass ``time_format = "original"`` (default = ``"reduced"``):

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
