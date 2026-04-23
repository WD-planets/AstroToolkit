"""
#########################
Working with Light Curves
#########################
Performing a Light Curve Query
==============================
Light curves can be fetched in a way that should be familiar from :doc:`previous tutorials <../getting_started/data_query>`. For this tutorial, ASAS-SN will be used:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits.gz")
asassn_query.show(show_types=True)
# sphinx_gallery_start_ignore
from ATK.queries.lightcurve._query_info import BAND_MAP

with open("../../auto_tutorials/lightcurves/supported_lightcurve_surveys.rst", "w") as f:
    for survey in BAND_MAP:
        f.write(f"        - {survey} ({', '.join(BAND_MAP[survey])})\n")
    f.write("\n")
# sphinx_gallery_end_ignore

# %%
# The returned :class:`~ATK.Models.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute is a list of :class:`Lightcurves <ATK.Models.Lightcurve>` (**one per photometric band per target, subject to data availability**).
#
# .. note::
# 
#    ATK supports queries to the following lightcurve surveys and bands:
# 
#    .. include:: supported_lightcurve_surveys.rst
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
#
# |
# |
#
# Plotting the Returned Data
# ==========================
# To plot a light curve and open it in the default browser:

# sphinx_gallery_start_ignore
asassn_query.plot()
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    For a refresher on plotting fundamentals, see :doc:`previous tutorials <../images/image_query>`.
#
# |
# |
#
# Splitting by Survey ID
# ======================
# **By default, light curves are returned as forced-photometry** (i.e. all detections in the requested radius are assumed to be uncontaminated by nearby sources). However, in some cases - e.g. in crowded regions or when using a large search radius - this assumption breaks down. To demonstrate this, we can change our target and increase the query radius to 1 arcminute:

import astropy.units as u

asassn_query = query("lightcurve", targets=587316166180416640, survey="asassn", radius = 60 * u.arcsec, path="example_lightcurve_2.fits.gz")
asassn_query.show(show_types=True)
# sphinx_gallery_start_ignore
asassn_query.plot()
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %% 
# It is clear that this is actually returning the photometry of at least two stars. In this case, it may be helpful to split the returned photometry by a per-survey observation/object ID. This can be done by passing ``split = True`` to :func:`~ATK.Tools.query`:

import astropy.units as u

asassn_query = query("lightcurve", targets=587316166180416640, survey="asassn",  radius = 60 * u.arcsec, split=True, path="example_lightcurve_3.fits.gz")
asassn_query.show(show_all=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    By default, :meth:`~ATK.Models.DataSet.show` truncates printing of the :attr:`~ATK.Models.DataSet.data` attribute if it contains a large number of containers. This can be disabled by passing ``show_all=True`` to :meth:`~ATK.Models.DataSet.show`, as above.
#
# |
#
# The now-split light curves can then be plotted in the same way as before - now producing a grid with three figures (one per unique source):

# sphinx_gallery_start_ignore
asassn_query.plot()
figure = format_plot(asassn_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
asassn_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. warning::
#
#    As ``split`` utilises a per-survey object ID to separate detections, its efficacy depends heavily on the specific implementation provided by each survey.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
