"""
###################
Working with Images
###################
To fetch an image, simply set the :func:`~ATK.Tools.query` ``kind`` to ``image`` and supply both a ``band`` and a ``size`` (rather than a ``radius``). An imaging ``survey`` must also be chosen - this tutorial will make use of Pan-STARRS (``panstarrs``):
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
import subprocess
# sphinx_gallery_end_ignore
from ATK import query

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, path="example_image_1.fits")
ps_query.show(show_types=True)
# sphinx_gallery_start_ignore
from ATK.queries.image._query_info import BAND_MAP
with open("../../auto_tutorials/images/supported_image_surveys.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    ATK supports queries to the following imaging surveys and bands:\n")
    for survey in BAND_MAP:
        f.write(f"        - {survey} - {', '.join(BAND_MAP[survey])}\n")
    f.write("\n")
# sphinx_gallery_end_ignore

# %%
# 
# As in any :func:`~ATK.Tools.query`, this returns a :class:`~ATK.Models.DataSet` with the :attr:`~ATK.Models.DataSet.data` attribute being a list of data containers - in this case :class:`~ATK.Models.Images`.
#
# | 
#
# .. note:: 
#
#    Since we have not supplied the units of ``size``, it has been assumed to be in arcsec. This can be changed in :doc:`the config <../configuration/config>`.
#
# .. include:: supported_image_surveys.rst
# 
# |
# |
#
# Plotting Data
# =============
# Unlike :class:`Records <ATK.Models.Record>`, :class:`~ATK.Models.Images` are plottable. The :meth:`~ATK.Models.DataSet.plot` method of the returned :class:`~ATK.Models.DataSet` can therefore be used to create a figure from the returned data:

ps_query.plot()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# The resulting figure will be saved to the :attr:`~ATK.Models.DataSet.figure` attribute of the :class:`~ATK.Models.DataSet`:

ps_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# The figure can be opened in the default browser by calling the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.open` method:

# sphinx_gallery_start_ignore
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
#
# .. note::
#
#    If a figure has not yet been generated when :meth:`~ATK.Models.DataSet.open` is called, :meth:`~ATK.Models.DataSet.plot` will be called automatically.
#
# |
# |
#
# Saving Figures to Local Files
# =============================
# To save a figure as a local HTML file, :meth:`~ATK.Models.DataSet.open` can be provided with a ``path``:
#
# .. code-block:: python
#
#    ps_query.open("example_image.html")
#
# .. warning::
#
#    If :meth:`~ATK.Models.DataSet.open` is not provided with a ``path``, the figure will be temporarily saved to the ``$HOME/.AstroToolkit/cached_figures`` directory. By default, cached figures that are over an hour old will be removed the next time :meth:`~ATK.Models.DataSet.open` is called.
#
# |
#
# Alternatively, the figure can be saved to local files without opening it by using the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.save` method:
#
# .. code-block:: python
#
#    ps_query.save("example_image.html")

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
