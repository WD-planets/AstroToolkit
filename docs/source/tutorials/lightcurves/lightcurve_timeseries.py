"""
###################
Timeseries Analysis
###################
Alongside the basic container methods shown in the :doc:`previous tutorial <lightcurve_manipulation>`, :class:`Lightcurves <ATK.Models.Lightcurve>` also support methods to facilitate :class:`Lomb-Scargle <astropy.timeseries.LombScargle>` timeseries analysis.

|

Generating Power Spectra
========================
To generate power spectra, we can use the :class:`~ATK.Models.Lightcurve.pspec` method to transform our lightcurves into :class:`~ATK.Models.Powspec` containers. By default, this process combines all bands for each target before computing a combined power spectrum.

Multi-band Power Spectra
------------------------
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits")
pspec_data = asassn_query.apply("pspec", min=0, max=60, samples=100000, inplace=False)
# sphinx_gallery_start_ignore
pspec_data.store("example_pspec.fits")
# sphinx_gallery_end_ignore
pspec_data.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# A :class:`~ATK.Models.Powspec` contains basic information about the light curve from which it was generated, along with :attr:`~ATK.Models.Powspec.frequency` and :attr:`~ATK.Models.Powspec.power` as numpy :class:`arrays <numpy.ndarray>`. The derived optimal frequency and the corresponding optimal period are stored in :attr:`~ATK.Models.Powspec.fopt` and :attr:`~ATK.Models.Powspec.popt`, respectively.
# 
# |
#
# We can then plot our power spectrum as with any other kind of data, by using :meth:`~ATK.Models.DataSet.plot` or :meth:`~ATK.Models.DataSet.open`:

# sphinx_gallery_start_ignore
pspec_data.plot()
figure = format_plot(pspec_data.figure, 1.5, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
pspec_data.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# Single-band Power Spectra
# -------------------------
# If we instead want to process each band individually, we can pass `multiband=False` to :meth:`~ATK.Models.Lightcurve.pspec`:

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits")
pspec_data = asassn_query.apply("pspec", min=0, max=60, samples=100000, multiband=False, inplace=False)
# sphinx_gallery_start_ignore
pspec_data.store("example_pspec.fits")
# sphinx_gallery_end_ignore
pspec_data.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# This generates one power spectrum per band:

# sphinx_gallery_start_ignore
pspec_data.plot()
figure = format_plot(pspec_data.figure, 1.5, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
pspec_data.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# Phase-Folding Light Curves
# ==========================
# Light curves can be phase folded onto a given period with the :meth:`~ATK.Models.Lightcurve.fold` method:

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits")
folded_data = asassn_query.apply("fold", min=0, max=60, samples=100000, multiband=False, inplace=False)
# sphinx_gallery_start_ignore
folded_data.store("example_folded_lc.fits")
# sphinx_gallery_end_ignore
folded_data.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# The :meth:`~ATK.Models.Lightcurve.fold` method simply replaces 

# sphinx_gallery_start_ignore
folded_data.plot()
figure = format_plot(folded_data.figure, 1.5, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open()
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
