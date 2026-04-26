"""
#####################
Multi-Target Plotting
#####################
Just like :func:`~ATK.Tools.query`, :meth:`~ATK.Models.DataSet.plot` can process data for multiple targets simultaneously. The following will perform a Pan-STARRS image query for four targets, before plotting and combining the returned data into a grid:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore

from ATK import query

targets = [2552928187080872832, 587316166180416640, 6050296829033196032, 4328198145165324800]
images = query("image", targets=targets, survey="panstarrs", size=120, band="g", path="multiple_images.fits.gz")
# sphinx_gallery_start_ignore
images.plot()
figure = format_plot(images.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
images.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %% 
# Plotting Arguments
# ==================
# Plotting arguments are applied to all plots in the grid:

# sphinx_gallery_start_ignore
images.plot(cmap="grey")
figure = format_plot(images.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
images.open(cmap="grey")
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
