"""
###################
Timeseries Analysis
###################
Alongside the basic data methods shown in the :doc:`previous tutorial <lightcurve_manipulation>`, :class:`Lightcurves <ATK.Models.Lightcurve>` also support :class:`Lomb-Scargle <astropy.timeseries.LombScargle>` timeseries analysis.

|

Generating Power Spectra
========================
The :class:`~ATK.Models.Lightcurve.pspec` method can be used to generate power spectra (i.e. :class:`~ATK.Models.Powspec` containers) from a set of :class:`Lightcurves <ATK.Models.Lightcurve>`. :class:`~ATK.Models.Lightcurve.pspec` must be provided with a minimum and maximum frequency, and a number of test frequencies in this range.

By default, this process combines all bands for each target before computing a combined power spectrum.

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

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits.gz")
pspec_data = asassn_query.apply("pspec", min=0, max=60, samples=100000, inplace=False)
# sphinx_gallery_start_ignore
pspec_data.store("example_pspec.fits.gz")
# sphinx_gallery_end_ignore
pspec_data.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# 
# |
#
# A :class:`~ATK.Models.Powspec` can be plotted like any other kind of data, with :meth:`~ATK.Models.DataSet.plot` or :meth:`~ATK.Models.DataSet.open`:

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
# |
# 
# Single-band Power Spectra
# -------------------------
# To instead process each band individually, pass ``multiband = False`` to :meth:`~ATK.Models.Lightcurve.pspec`:

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits.gz")
pspec_data = asassn_query.apply("pspec", min=0, max=60, samples=100000, multiband=False, inplace=False)
# sphinx_gallery_start_ignore
pspec_data.store("example_pspec.fits.gz")
# sphinx_gallery_end_ignore
pspec_data.show(show_types=True)
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
# 
# |
# |
#
# Phase-Folding Light Curves *
# ============================
# :class:`Lightcurves <ATK.Models.Lightcurve>` can be phase folded on a given period with the :meth:`~ATK.Models.Lightcurve.fold` method:

asassn_query = query("lightcurve", targets=6050296829033196032, survey="asassn", path="example_lightcurve.fits.gz")
folded_data = asassn_query.apply("fold", min=0, max=60, samples=100000, multiband=False, inplace=False)
# sphinx_gallery_start_ignore
folded_data.store("example_folded_lc.fits.gz")
# sphinx_gallery_end_ignore
folded_data.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%

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
