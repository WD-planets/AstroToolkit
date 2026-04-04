"""
#################
SED Customisation
#################
Overlaying Spectra
==================
SEDs can be overlayed with spectra by passing a :class:`~ATK.Models.DataSet` containing :class:`Spectra <ATK.Models.Spectrum>` to :meth:`~ATK.Models.DataSet.plot` as an ``overlay``:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

spec_query = query("spectrum", targets=587316166180416640, survey="sdss", path="example_spectrum.fits.gz")
sed_query = query("sed", targets=587316166180416640, path="example_sed.fits.gz")

# sphinx_gallery_start_ignore
sed_query.plot(overlay=spec_query)
figure = format_plot(sed_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sed_query.open(overlay=spec_query)
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    :class:`Spectra <ATK.Models.Spectrum>` and their corresponding :class:`SEDs <ATK.Models.SED>` are automatically matched by :class:`~ATK.Models.Target` (i.e. by input search position or Gaia source ID) before they are overlayed. In the case of multiple matching :class:`Spectra <ATK.Models.Spectrum>`, each one is overlayed over a copy of the :class:`~ATK.Models.SED`.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
