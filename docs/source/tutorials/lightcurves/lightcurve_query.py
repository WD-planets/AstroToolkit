"""
#########################
Working with Light Curves
#########################

We can fetch light curves in a way that is now hopefully familiar:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits")
asassn_query.show(show_types=True)
# sphinx_gallery_start_ignore
from ATK.queries.lightcurve._query_info import BAND_MAP

with open("../../auto_tutorials/lightcurves/supported_lightcurve_surveys.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    ATK supports queries to the following lightcurve surveys and bands:\n")
    for survey in BAND_MAP:
        f.write(f"        - {survey} - {', '.join(BAND_MAP[survey])}\n")
    f.write("\n")
# sphinx_gallery_end_ignore

# %%
# This returns a :class:`~ATK.Models.DataSet` with the :attr:`~ATK.Models.DataSet.data` attribute being a list of returned :class:`~ATK.Models.Lightcurve` containers (**one per photometric band per target**). Each :class:`~ATK.Models.Lightcurve` contains information about the search, including the ``separation`` between the light curve and the position of the search, along with the actual photometric data as numpy :class:`arrays <numpy.ndarray>`.
#
# |
#
# .. note::
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
# 
# .. include:: supported_lightcurve_surveys.rst
#
# |
# |
#
# Plotting Data
# =============
# To plot a basic light curve and open it in the default browser:

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
# By default, light curves are returned as forced-photometry (i.e. all detections in the requested radius are assumed to be uncontaminated by nearby sources). However, in crowded regions this may not be true. To demonstrate this, we can change our target and increase the query radius to 1 arcminute:

import astropy.units as u

asassn_query = query("lightcurve", targets=587316166180416640, survey="asassn", path="example_lightcurve.fits", radius = 60 * u.arcsec, split=False)
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
# It is clear that we are actually getting the light curves of at least two systems. In this case, it may help to split the returned photometry by a per-survey observation/object ID. This can be done by passing ``split = True`` to :func:`~ATK.Tools.query`:

import astropy.units as u

asassn_query = query("lightcurve", targets=587316166180416640, survey="asassn", path="example_lightcurve.fits", radius = 60 * u.arcsec, split=True)
asassn_query.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    By default, :meth:`~ATK.Models.DataSet.show` will truncate printing when a :class:`~ATK.Model.DataSet` contains many data containers. This can be disabled by passing ``show_all=True`` to :meth:`~ATK.Models.DataSet.show`, as above.
#
# This splits the returned :class:`Lightcurves <ATK.Models.Lightcurve>` into whatever the survey considers to be distinct detections. In this case, we can see that our original light curve has been split into three (hence the three unique ``separations`` of our returned :class:`Lightcurves <ATK.Models.Lightcurve>`).
#
# We can then :meth:`~ATK.Models.DataSet.plot` our newly split :class:`Lightcurves <ATK.Models.Lightcurve>` as usual:

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
#
# Multi-Target Plotting
# =====================
# If a :class:`~ATK.Models.DataSet` contains data for multiple targets, plotting will automatically sort thjese

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
