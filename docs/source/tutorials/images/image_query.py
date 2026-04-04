"""
###################
Working with Images
###################
Performing an Image Query
=========================
To fetch an image, set the :func:`~ATK.Tools.query` ``kind`` to ``"image"`` and supply both a ``band`` and a ``size`` (rather than a ``radius``). A ``survey`` must also be chosen - this tutorial will use Pan-STARRS:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, path="example_image_1.fits.gz")
ps_query.show(show_types=True)
# sphinx_gallery_start_ignore
from ATK.queries.image._query_info import BAND_MAP
with open("../../auto_tutorials/images/supported_image_surveys.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    ATK supports queries to the following imaging surveys and bands:\n")
    for survey in BAND_MAP:
        f.write(f"        - {survey} ({', '.join(BAND_MAP[survey])})\n")
    f.write("\n")
    f.write("    If the desired survey is not listed above, see :doc:`here <../extension/external_data>` for a tutorial on utilising external data.")
# sphinx_gallery_end_ignore

# %%
# 
# As in any :func:`~ATK.Tools.query`, this returns a :class:`~ATK.Models.DataSet`. In this case, the :class:`~ATK.Model.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute is a list of :class:`~ATK.Models.Images` (**one per target, subject to data availability**).
#
# | 
#
# .. include:: supported_image_surveys.rst
#
# |
# |
#
# Plotting the Returned Data
# ==========================
# :class:`Images <ATK.Models.Image>` are plottable. The returned :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.plot` method can therefore be used to create a Bokeh :class:`~bokeh.plotting.figure`:

ps_query.plot()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# This saves a  figure to the :attr:`~ATK.Models.DataSet.figure` attribute of the :class:`~ATK.Models.DataSet`:

ps_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
#
# |
# |
#
# Viewing a Figure
# ================
# A :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.open` method can be used to open the stored :class:`~bokeh.plotting.figure` in the default browser:

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
# **Since a Gaia source ID was used for targeting, proper motion correction has been used to automatically centre the image on the target star.**
#
# .. note::
#
#    If a figure has not yet been generated when :meth:`~ATK.Models.DataSet.open` is called, :meth:`~ATK.Models.DataSet.plot` is called automatically.
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
#    If :meth:`~ATK.Models.DataSet.open` is not provided with a ``path``, the figure will be temporarily saved to the ``~/.AstroToolkit/cached_figures`` directory (i.e. inside the home directory). By default, cached figures that are over an hour old will be removed the next time :meth:`~ATK.Models.DataSet.open` is called.
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
