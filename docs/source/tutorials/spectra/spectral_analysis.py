"""
#################
Analysing Spectra
#################

Velocity Spectra
================
:class:`Spectra <ATK.Models.DataSet>` support one **data method**, :meth:`~ATK.Models.Spectrum.vspec`, which converts a spectrum into a velocity spectrum relative to a given reference wavelength, ``wav_ref``:
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
# :class:`Spectra <ATK.Models.DataSet>` also support two **plot methods**. These do not modify the :class:`~ATK.Models.Spectrum` containers themselves, but instead add additional elements to plots.
#
# .. note:: 
# 
#    If a plot has not been generated at the time a **plot method** is called, one will be automatically generated with default parameters.
# 
# |
#
# Spectral features can be detected and fitted with :meth:`~ATK.Models.Spectrum.fit`. To aid in detecting features while reducing false positives, :meth:`~ATK.Models.Spectrum.fit` supports the following arguments:
# 
# - ``prominence`` sets the minimum **prominence** of spectral features (see :func:`here <scipy.signal.find_peaks>` for details, default = ``2.0``)
# 
# - ``smooth`` sets the level to which the spectrum is **smoothed** before searching for peaks (default = ``3.0``)
# 
# - ``snr`` sets the minimum **signal-to-noise ratio** of spectral features (default = ``3.0``)
# 
# |

sdss_query.apply("fit", smooth=5.0)
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
# :class:`Spectra <ATK.Models.Spectrum>` also support multi-component radial velocity fitting via :meth:`~ATK.Models.Spectrum.rv_fit`. As peaks must first be detected, :meth:`~ATK.Models.Spectrum.rv_fit` accepts the same parameters as :meth:`~ATK.Models.Spectrum.fit`:

sdss_query.apply("rv_fit", smooth=5.0)
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
