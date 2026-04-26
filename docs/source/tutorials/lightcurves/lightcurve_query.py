"""
#########################
Working with Light Curves
#########################
Performing a Light Curve Query
==============================
Light curves can be fetched in the same way as in :doc:`previous tutorials <../getting_started/data_query>`. For this tutorial, an ASAS-SN query will be performed:
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
        f.write(f"- {survey} ({', '.join(BAND_MAP[survey])})\n")
    f.write("\n")
# sphinx_gallery_end_ignore

# %%
# |
#
# The returned :class:`~ATK.Models.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute is a list of :class:`~ATK.Models.Lightcurve` objects (**one per photometric band per target, subject to data availability**).
# 
# |
#
# .. note::
# 
#    ATK supports queries to the following lightcurve surveys and bands:
# 
#    .. include:: supported_lightcurve_surveys.rst
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
# 
#    If the desired survey is not listed above, see :doc:`here <../extension/external_data>` for a tutorial on utilising external data.
#
# |
# |
#  
# Plotting the Returned Data
# ==========================
# To plot a :class:`~ATK.Models.Lightcurve` and open it in the default browser:

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
# **By default, light curves are returned as forced-photometry** (i.e. all photometric detections in the requested radius are returned as a single :class:`~ATK.Models.Lightcurve` per band). However, in some cases - e.g. in crowded regions or when using a large search radius - this assumption breaks down:

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
# |
#
# **It is clear that the returned** :class:`~ATK.Models.Lightcurve` **objects contain photometry of at least two stars.** In this case, it may be helpful to split the returned photometry by a per-survey observation/object ID. This can be done by passing ``split=True`` to :func:`~ATK.Tools.query`:

import astropy.units as u

asassn_query = query("lightcurve", targets=587316166180416640, survey="asassn",  radius = 60 * u.arcsec, split=True, path="example_lightcurve_3.fits.gz")
asassn_query.show(show_all=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# |
# 
# .. note::
#
#    By default, :meth:`~ATK.Models.DataSet.show` truncates printing of the :attr:`~ATK.Models.DataSet.data` attribute if it contains a large number of containers. This can be disabled by passing ``show_all=True`` to :meth:`~ATK.Models.DataSet.show`.
#
# |
#
# The split light curves can then be plotted in the same way as before - now producing a grid with three figures (one per unique source):

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
#
# |
# |
#
# Data Quality Filtering
# ======================
# By default, basic quality filtering is enabled for the following surveys:
#
#     **ATLAS**: Filtering is performed as advised `here <https://fallingstar-data.com/forcedphot/faq/>`_.
#
#     **TESS**: ``Flux significance > 3``, ``quality flag = 0``.
#
# |
#
# To disable all **unrequired** filtering, pass ``filter=False`` to :func:`~ATK.Tools.query` when retrieving data from one the above surveys.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
