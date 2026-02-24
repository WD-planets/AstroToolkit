"""
############################
Fetching and Plotting Images
############################
To fetch an image, we simply set the :func:`~ATK.Tools.query` ``kind`` to ``image`` and supply both a ``band`` and a ``size`` (rather than a ``radius``). We must also choose an imaging survey - in this tutorial we will use Pan-STARRS (``panstarrs``):
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
ps_query.show()
# sphinx_gallery_start_ignore
from ATK.queries.image._query_info import BAND_MAP
with open("../../auto_tutorials/images/supported_image_surveys.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    ATK supports queries to the following imaging surveys:\n")
    for survey in BAND_MAP:
        f.write(f"        - {survey} - {', '.join(BAND_MAP[survey])}\n")
    f.write("\n")
# sphinx_gallery_end_ignore

# %%
# 
# As with all queries, this returns a :class:`~ATK.Models.DataSet` with the ``data`` attribute being a list of data containers - in this case :class:`~ATK.Models.Images`.
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
# Unlike :class:`VizierEntries <ATK.Models.VizierEntry>`, :class:`~ATK.Models.Images` are plottable. This means that we can use the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.plot` method to create a figure from the returned data:

ps_query.plot()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# The resulting figure will be saved to the :attr:`~ATK.Models.DataSet.figure` attribute of the :class:`~ATK.Models.DataSet`:

ps_query.show()

# %%
# We can then open the figure in the default browser by calling the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.open` method:

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
# 
# If we want to also save this figure as a local HTML page, we can provide :meth:`~ATK.Models.DataSet.open` with a ``path``:
# 
# .. code-block:: python
# 
#    ps_query.open("example_image.html")
#
# .. note::
#    
#    If :meth:`~ATK.Models.DataSet.open` is not provided with a ``path``, the figure will be temporarily saved to the ``$HOME/.AstroToolkit/cached_figures`` directory. By default, cached figures that are over an hour old will be removed the next time :meth:`~ATK.Models.DataSet.open` is called.
# 
# |
# 
# Alternatively, we can save the figure to local files without opening it by using the :meth:`~ATK.Models.DataSet.save` method:
# 
# .. code-block:: python
# 
#    ps_query.save("example_image.html")
# 
# |
# |
# 
# Detection Overlays
# ==================
# For catalogues that are supported out-of-the-box (see :doc:`the first tutorial <../getting_started/data_query>`), detection overlays are already configured. We can therefore immediately request a GALEX overlay, which is done by passing the ``overlays`` parameter to :func:`~ATK.Tools.query`: 

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, overlays=["galex"], path="example_image_2.fits")
# sphinx_gallery_start_ignore
ps_query.plot()
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
#    Detection markers are corrected for proper motion by utilising the catalogue's epoch (as set in the epoch file). Non-Gaia detections are corrected where possible by "piggybacking" them with corresponding Gaia detections within a radius. This radius can be set in :doc:`the config <../configuration/config>`. 
#
# | 
#
# Defining a Detection Overlay
# ----------------------------
# For user-added catalogues (e.g. our AllWISE alias from the :doc:`previous tutorials <../getting_started/data_query_2>`), we need to provide ATK with the necessary information to generate a detection overlay. This is done via the overlay file, which can be accessed and edited in the same way as the alias and epoch files:
#
# .. code-block:: console
#
#    $ ATKoverlay show

# %% 
# from the command line, or:

from ATK.Config import OVERLAY_CONFIG

# sphinx_gallery_start_ignore
OVERLAY_CONFIG.reset()
# sphinx_gallery_end_ignore
OVERLAY_CONFIG.show()

# %%
# from inside a script. We can then use:
#
# .. code-block:: console
#
#    $ ATKoverlay open
#
# or: 
#
# .. code-block:: python
#
#    from ATK.Config import OVERLAY_CONFIG
#
#    OVERLAY_CONFIG.open()

# %%
# to add a ``photometric`` overlay definition:
#
# .. code-block::
#
#    allwise
#        mags: [W1mag, W2mag, W3mag, W4mag]
#        errors: [e_W1mag, e_W2mag, e_W3mag, e_W4mag]
#        lon_column = RAJ2000
#        lat_column = DEJ2000
#        frame = icrs
#
# .. note::
# 
#    If we only care about the positions of any detections, we can instead add a ``positional`` overlay definition without supplying ``mags`` or ``errors``.
#
# |
# 
# We can now request a GALEX and AllWISE detection overlay:

# sphinx_gallery_start_ignore
subprocess.run(["ATKoverlay", "set", "photometric", "allwise", "--mags", "W1mag", "W2mag", "W3mag", "W4mag", "--errors", "e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag", "--lon", "RAJ2000", "--lat", "DEJ2000", "--frame", "icrs"])
# sphinx_gallery_end_ignore

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, overlays=["galex", "allwise"], path="example_image_3.fits")
# sphinx_gallery_start_ignore
ps_query.plot()
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
# |
#
# Configuring an Overlay 
# ----------------------
# By default, the detection overlay uses whichever band is listed first in the catalogue's overlay definition. If we wish to overlay a different band, or perhaps overlay multiple bands from the same catalogue simultaneously, we can instead pass overlays as a ``dict``:

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, overlays={"galex":["NUVmag","FUVmag"], "allwise":["W4mag"]}, path="example_image_4.fits")
# sphinx_gallery_start_ignore
ps_query.plot()
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
ps_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# .. note::
#
#    Overlayed detection markers are scaled in size by magnitude (larger = brighter). Markers labelled ``galex detection`` in the above image are detections that had an invalid (i.e. ``NaN``) magnitude, and hence could not be scaled.

# %%
#
# |
# |
# |
#
# Download this Tutorial
# ----------------------
