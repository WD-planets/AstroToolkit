"""
#################
Analysing Spectra
#################

Velocity Spectra
================
:class:`~ATK.Models.DataSet` objects support one **data method**, :meth:`~ATK.Models.Spectrum.vspec`, which converts wavelength into velocity relative to a given reference wavelength, ``wav_ref``:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

sdss_query = query("spectrum", targets=587316166180416640, survey="sdss", path="example_spectrum.fits.gz")
sdss_query.apply("vspec", wav_ref = 6562.7097, inplace=False)
# sphinx_gallery_start_ignore
sdss_query.plot()
figure = format_plot(sdss_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sdss_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
# |
# 
# Fitting Spectral Features
# =========================
# :class:`~ATK.Models.Spectrum` objects also support two further plotting options for performing spectral analysis.
#
# |
#
# **Spectral features can be detected and fitted** by passing ``fit=True`` to :meth:`~ATK.Models.DataSet.plot`. To aid in detecting features while reducing false positives, the following additional arguments are also supported when using ``fit=True``:
# 
# - ``prominence`` sets the minimum **prominence** of spectral features (see :func:`here <scipy.signal.find_peaks>` for details, default = ``2.0``)
# 
# - ``smooth`` sets the level to which the spectrum is **smoothed** before searching for peaks (default = ``3.0``)
# 
# - ``snr`` sets the minimum **signal-to-noise ratio** of spectral features (default = ``3.0``)
# 
# |

sdss_query.plot(fit=True, smooth=5.0)
# sphinx_gallery_start_ignore
figure = format_plot(sdss_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sdss_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
# |
#
# Detecting Radial Velocities
# ===========================
# :class:`~ATK.Models.Spectrum` plotting also supports **multi-component radial velocity fitting** by passing ``rv_fit=True``. As peaks must first be detected, **the same additional parameters are accepted as when using** ``fit=True``:

sdss_query.plot(rv_fit=True, smooth=5.0)
# sphinx_gallery_start_ignore
figure = format_plot(sdss_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sdss_query.open()
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
