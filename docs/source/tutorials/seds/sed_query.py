"""
#################
Working with SEDs
#################
To perform an SED query and plot the result:
"""

# sphinx_gallery_start_ignore
from ATK.queries.sed._query_info import SURVEY_MAP

with open("../../auto_tutorials/seds/supported_sed_surveys.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    ATK generates SEDs using photometry from the following surveys:\n")
    for survey in SURVEY_MAP:
        f.write(f"        - {survey}\n")
    f.write("\n")
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

sed_query = query("sed", targets=587316166180416640, path="example_sed.fits")
sed_query.show(show_types=True)

# sphinx_gallery_start_ignore
sed_query.plot()
figure = format_plot(sed_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sed_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# This returns a :class:`~ATK.Models.DataSet` with the :attr:`~ATK.Models.DataSet.data` attribute being a list of returned :class:`~ATK.Models.SED` containers (one per target). Each :class:`~ATK.Models.SED` contains all retrieved photometry within the search radius.
#
# |
#
# .. note::
#
#    For a refresher on :func:`~ATK.Tools.query` fundamentals, see :doc:`previous tutorials <../getting_started/data_query>`.
#
# .. include:: supported_sed_surveys.rst
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
