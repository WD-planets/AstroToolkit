"""
Defining Custom Image Overlays
==============================
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore

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

import subprocess

from ATK.Config import SURVEY_CONFIG
from ATK import query

# sphinx_gallery_start_ignore
SURVEY_CONFIG.reset()
# sphinx_gallery_end_ignore
SURVEY_CONFIG.show()

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
        "ATKsurvey",
        "set",
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

ps_query = query(
    "image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, overlays=["galex", "allwise"], path="example_image_3.fits.gz"
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
#    If photometry is not needed/is not present, a ``positional`` overlay definition can instead be defined without supplying ``mags`` or ``errors``.
#
# |
