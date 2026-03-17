"""
#################
Working with HRDs
#################
To perform a HRD query and plot the result:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

hrd_query = query("hrd", targets=[587316166180416640], path="example_hrd.fits")
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
# This returns a :class:`~ATK.Models.DataSet` with the :attr:`~ATK.Models.DataSet.data` attribute being a list of returned :class:`~ATK.Models.HRD` containers. <MULTIPLE CONTAINERS>
# |
#
# .. note::
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
