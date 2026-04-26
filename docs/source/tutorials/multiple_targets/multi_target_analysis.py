"""
#####################
Multi-Target Analysis
#####################
As with plotting, analysis can be performed in exactly the same way as with a single target - :meth:`~ATK.Models.DataSet.apply` will apply all **data methods** across all containers in the :class:`~ATK.Models.DataSet`:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

targets = [5346631922949364864, 6050296829033196032]
asassn_query = query("lightcurve", targets=targets, survey="asassn", path="multiple_lightcurves.fits.gz")
folded_data = asassn_query.apply("fold", fmin=0, fmax=60, samples=100000, inplace=False)
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
