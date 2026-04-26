"""
#################
Working with HRDs
#################
Performing a HRD Query
======================
To perform a Gaia HRD query and plot the result:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

hrd_query = query("hrd", targets=587316166180416640, path="example_hrd.fits.gz")
hrd_query.show(show_types=True)

# sphinx_gallery_start_ignore
hrd_query.plot()
figure = format_plot(hrd_query.figure, 2, 2)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
hrd_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# The returned :class:`~ATK.Models.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute is a list of :class:`~ATK.Models.HRD` objects (**one per target, subject to data availability**).
#
# |
#
# .. note::
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`here <../getting_started/data_query>`. For a refresher on plotting fundamentals, see :doc:`here <../images/image_query>`.
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
