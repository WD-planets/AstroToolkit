"""
#############################
Working with Multiple Surveys
#############################
Combining Data Sets
===================
The :meth:`~ATK.Models.DataSet.merge` method of a :class:`~ATK.Models.DataSet` can be used to combine it with another, **as long as they store data of the same type** (for a tutorial on plotting multiple types of data together, see :doc:`here <../datapages/datapage_plotting>`). Along with making the recombination of :meth:`~ATK.Models.Dataset.split` :class:`DataSets <ATK.Models.DataSet>` possible, this also allows for data from more than one survey to be utilised simultaneously:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

target = 587316166180416640

asassn_data = query("lightcurve", targets=target, survey="asassn", path="asassn_lightcurve.fits.gz")
gaia_data = query("lightcurve", targets=target, survey="gaia", path="gaia_lightcurve.fits.gz")

merged_data = asassn_data.merge(gaia_data)
merged_data.show()
# sphinx_gallery_start_ignore
merged_data.plot()
figure = format_plot(merged_data.figure, 3, 1.5, False)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
merged_data.open()
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
