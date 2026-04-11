"""
#################
HRD Customisation
#################
Combining Multiple HRDs
=======================
By default, plotting a :class:`~ATK.Models.DataSet` with multiple :class:`HRDs <ATK.Models.ATK>` generates multiple plots. To instead combine all sources into a single plot, pass ``split = False`` to :meth:`~ATK.Models.DataSet.plot` or :meth:`~ATK.Models.DataSet.open`:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

hrd_query = query("hrd", targets=[587316166180416640, 6050296829033196032], path="example_hrd_2.fits.gz")

# sphinx_gallery_start_ignore
hrd_query.plot(split=False)
figure = format_plot(hrd_query.figure, 2, 2)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
hrd_query.open(split=False)
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# 
# .. note:: 
# 
#    The above is an example of a multi-target query. For more information, see the :doc:`next section <../multiple_targets/multi_target_query>`.
#
# |
#
# Reducing the Background
# =======================
# Due to the number of background data points, plotted :class:`HRDs <ATK.Models.HRD>` can be expensive in terms of browser performance. This can be avoided by passing ``background`` as a ``float`` between ``0`` and ``1`` (default = ``1``):

hrd_query = query("hrd", targets=587316166180416640, path="example_hrd.fits.gz")

# sphinx_gallery_start_ignore
hrd_query.plot(background=0.25)
figure = format_plot(hrd_query.figure, 2, 2)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
hrd_query.open(background=0.25)
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
