"""
####################
Working with Spectra
####################
Performing a Spectrum Query
===========================
To perform a spectrum query and plot the result:
"""

# sphinx_gallery_start_ignore
from ATK.queries.spectrum._query_info import SURVEY_MAP

with open("../../auto_tutorials/spectra/supported_spectrum_surveys.rst", "w") as f:
    for survey in SURVEY_MAP:
        f.write(f"    - {survey}\n")
    f.write("\n")
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

sdss_query = query("spectrum", targets=587316166180416640, survey="sdss", path="example_spectrum.fits.gz")
sdss_query.show(show_types=True)
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
# The returned :class:`~ATK.Models.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute is a list of :class:`~ATK.Models.Spectrum` containers (**one or multiple per target, subject to data availability**). If multiple :class:`Spectra <ATK.Models.Spectrum>` are returned for a single target, these will be arranged in order of decreasing exposure.
#
# |
#
# .. note::
#
#    ATK supports queries to the following spectra surveys:
#
#    .. include:: supported_spectrum_surveys.rst
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
#
#    If the desired survey is not listed above, see :doc:`here <../extension/external_data>` for a tutorial on utilising external data.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
