"""
###################
Timeseries Analysis
###################
Alongside the basic **data methods** shown in the :doc:`previous tutorial <lightcurve_manipulation>`, :class:`~ATK.Models.Lightcurve` objects also support Lomb-Scargle timeseries analysis.

|

Generating Power Spectra
========================
The :meth:`~ATK.Models.Lightcurve.pspec` method can be used to generate power spectra from one or more :class:`Lightcurves <ATK.Models.Lightcurve>`. A minimum and maximum frequency must be provided, along with a number of test frequencies in this range. **If no units are given, values are assumed to be in** :math:`\mathbf{days}^{-1}`.

This process transforms the :class:`~ATK.Models.DataSet`'s :class:`~ATK.Models.Lightcurve` objects into one or more :class:`~ATK.Models.Powspec` objects.

Multi-band Power Spectra
------------------------
By default, :meth:`~ATK.Models.Lightcurve.pspec` combines all bands for each target before computing a combined power spectrum.
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=5346631922949364864, survey="asassn", path="example_lightcurve.fits.gz")
pspec_data = asassn_query.apply("pspec", fmin=0, fmax=10, samples=10000, inplace=False)
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
# To instead process each band individually, pass ``multiband = False`` to :meth:`~ATK.Models.Lightcurve.pspec`. This treats each band as being entirely separate:

asassn_query = query("lightcurve", targets=5346631922949364864, survey="asassn", path="example_lightcurve.fits.gz")
pspec_data = asassn_query.apply("pspec", fmin=0, fmax=10, samples=10000, multiband=False, inplace=False)
# sphinx_gallery_start_ignore
pspec_data.store("example_pspec.fits.gz")
# sphinx_gallery_end_ignore
pspec_data.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# |
#
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
# Phase-Folding Light Curves
# ==========================
# :class:`~ATK.Models.Lightcurve` objects can be phase folded with :meth:`~ATK.Models.Lightcurve.fold`. **By default,** :meth:`~ATK.Models.Lightcurve.fold` **first generates a power spectrum from which to derive an optimal frequency.** :meth:`~ATK.Models.Lightcurve.fold` therefore accepts all parameters that can be passed to :meth:`~ATK.Models.Lightcurve.pspec`.
#
# A multiband phase-folded light curve can threfore be generated and plotted as follows:

asassn_query = query("lightcurve", targets=5346631922949364864, survey="asassn", path="example_lightcurve.fits.gz")
folded_data = asassn_query.apply("fold", fmin=0, fmax=10, samples=10000, inplace=False)
# sphinx_gallery_start_ignore
folded_data.store("example_folded_lc.fits.gz")
# sphinx_gallery_end_ignore
folded_data.show(show_types=True)
# sphinx_gallery_start_ignore
folded_data.plot()
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
# 
# Frequency Optimisation
# ----------------------
# Despite the underlying power spectrum being the same, the photometry in the above example has been folded on a different peak frequency - in this case the first subharmonic (i.e. ``0.5 * fopt``). 
# 
# **By default, a phase-dispersion metric is calculated for a set of harmonics either side of the peak frequency** and the best is chosen. This can help to preserve real periodic structure, especially if the modulation is asymmetric as seen here. 
# 
# This behaviour can be disabled by passing ``optimise=False`` to :meth:`~ATK.Models.DataSet.apply` - in this case losing the true orbital frequency in favour of the peak frequency in the power spectrum:

folded_data = asassn_query.apply("fold", fmin=0, fmax=10, samples=10000, optimise=False, inplace=False)
# sphinx_gallery_start_ignore
folded_data.store("example_folded_lc.fits.gz")
# sphinx_gallery_end_ignore
folded_data.show(show_types=True)
# sphinx_gallery_start_ignore
folded_data.plot()
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Frequency Overriding
# --------------------
# To avoid the derivation of an optimal frequency entirely, a frequency can be passed via the ``freq`` parameter:

folded_data = asassn_query.apply("fold", freq=2.773227732277323, inplace=False)
# sphinx_gallery_start_ignore
folded_data.store("example_folded_lc.fits.gz")
# sphinx_gallery_end_ignore
folded_data.show(show_types=True)
# sphinx_gallery_start_ignore
folded_data.plot()
figure = format_plot(folded_data.figure, 3, 1.5)
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
