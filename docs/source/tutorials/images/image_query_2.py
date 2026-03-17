"""
##################
Detection Overlays
##################
To simplify the matching of positional data, :class:`Images <ATK.Models.Image>` support detection overlays. The following catalogues are already configured for this by default:

.. include:: ../getting_started/supported_aliases.rst

An overlay can be requested by passing ``overlays`` to :func:`~ATK.Tools.query`:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
import subprocess
# sphinx_gallery_end_ignore
from ATK import query

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
#    Detection markers are corrected for proper motion by utilising the catalogue's epoch from the epoch file. Non-Gaia detections are corrected where possible by "piggybacking" them with nearby Gaia detections.
#
# |
# |
#
# Defining a Detection Overlay
# ============================
# For user-added catalogue aliases (e.g. the one added in a :doc:`previous tutorial <data_query_2>`), ATK needs to be supplied with the necessary information to generate a detection overlay.  
#
# |
# 
# Adding an Overlay Definition
# ----------------------------
# A new overlay definition can be added via the overlay file, which can be accessed and edited in the same way as the alias and epoch files seen in a :doc:`previous tutorial <data_query_2>`:
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
# from inside a script. By opening the file with:
#
# .. code-block:: console
#
#    $ ATKoverlay open
#
# .. code-block:: python
#
#    from ATK.Config import OVERLAY_CONFIG
#
#    OVERLAY_CONFIG.open()

# %%
# a ``photometric`` overlay definition can be added:
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
# |
#
# Testing the New Overlay 
# ------------------------
# An AllWISE detection overlay can now be requested alongside our GALEX overlay:

# sphinx_gallery_start_ignore
subprocess.run(
    [
        "ATKoverlay",
        "set",
        "photometric",
        "allwise",
        "--mags",
        "W1mag",
        "W2mag",
        "W3mag",
        "W4mag",
        "--errors",
        "e_W1mag",
        "e_W2mag",
        "e_W3mag",
        "e_W4mag",
        "--lon",
        "RAJ2000",
        "--lat",
        "DEJ2000",
        "--frame",
        "icrs",
    ]
)
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
# .. note::
#
#    If photometry is not needed/is not present, a ``positional`` overlay definition can instead be defined without supplying ``mags`` or ``errors``.
#
# |
#
# Further Overlay Configuration
# -----------------------------
# By default, a detection overlay uses whichever band is listed first in the catalogue's overlay definition. If a different band is required, or perhaps multiple bands from the same catalogue need to be overlayed simultaneously, ``overlays`` can instead be passed as a ``dict``:

ps_query = query(
    "image",
    targets=2552928187080872832,
    survey="panstarrs",
    band="g",
    size=120,
    overlays={"galex": ["NUVmag", "FUVmag"], "allwise": ["W4mag"]},
    path="example_image_4.fits",
)
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
# .. rubric:: Download this Tutorial
