"""
###################
Image Customisation
###################
Similarly to :func:`~ATK.Tools.query`, the :meth:`~ATK.Models.DataSet.plot` method of a :class:`~ATK.Models.DataSet` accepts various additional arguments depending on the kind of data that is being plotted.

When plotting an :class:`~ATK.Models.Image`, two such arguments can be utilised: ``cmap`` and ``relative_axes``.

|

.. note::

   Calling :meth:`~ATK.Models.DataSet.open` automatically calls :meth:`~ATK.Models.DataSet.plot` if a figure has not already been generated. In this case, :meth:`~ATK.Models.DataSet.open` passes any additional keyword arguments to :meth:`~ATK.Models.DataSet.plot`.
   A new plot will also be generated if the keyword arguments passed to :meth:`~ATK.Models.DataSet.open` do not match those of the stored figure.

|

Changing the Colour Map
=======================
The ``cmap`` parameter sets the colour map of the plotted image:

1. ``viridis`` (default):
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore

from ATK import query

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, path="example_image_1.fits.gz")
# sphinx_gallery_start_ignore
ps_query.plot(cmap="viridis")
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open(cmap="viridis")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# 
# |
# 
# 2. ``grey``:

# sphinx_gallery_start_ignore
ps_query.plot(cmap="grey")
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open(cmap="grey")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
#
# |
#
# 3. ``false_colour``, which maps the wavelength of the image filter into a single RGB colour:

# sphinx_gallery_start_ignore
ps_query.plot(cmap="false_colour")
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open(cmap="false_colour")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
#
# |
# |
#
# Configuring Axes Coordinates
# ============================
# So far, all images have had coordinate axes that are defined relative to the image centre, i.e. ``relative_axes = True`` (default). To instead use absolute coordinates on the sky, `pass ``relative_axes = False``:

# sphinx_gallery_start_ignore
ps_query.plot(relative_axes=False)
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open(relative_axes=False)
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
