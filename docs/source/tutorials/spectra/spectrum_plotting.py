"""
######################
Spectrum Customisation
######################
Overlaying SEDs
===============
Spectra can be overlayed with a spectral energy distribution (see the :doc:`next tutorial <../seds/sed_query>`) by passing a :class:`~ATK.Models.DataSet` containing :class:`~ATK.Models.SED` objects to :meth:`~ATK.Models.DataSet.plot` or :meth:`~ATK.Models.DataSet.open`:
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
spec_query.plot(overlay=sed_query)
figure = format_plot(spec_query.figure, 3, 2)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
spec_query.open(overlay=sed_query)
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# .. note::
#
#    :class:`~ATK.Models.Spectrum` objects and their corresponding :class:`~ATK.Models.SED` overlays are automatically matched before being overlayed.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
