"""
##################
Detection Overlays
##################
To aid the matching of positional data, image queries support detection overlays.

.. note::

   By default, detection overlays are configured for the following Vizier catalogues:

   .. include:: ../getting_started/supported_aliases.rst

   A tutorial on defining detection overlays for any Vizier catalogue can be found :doc:`here <../extension/overlays>`.

|
|

Requesting an Overlay
=====================
An overlay can be requested by passing ``overlays`` to :func:`~ATK.Tools.query` with a list of surveys:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

ps_query = query("image", targets=2552928187080872832, survey="panstarrs", band="g", size=120, overlays=["galex"], path="example_image_2.fits")
# sphinx_gallery_start_ignore
ps_query.plot()
figure = format_plot(ps_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore'/home/ethan/Documents/atk_1.8/ATK_new/docs/build/html/tutorials/seds/custom_index.html' 
ps_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# 
# As with all figure elements in ATK, overlayed detections can be hidden by clicking them in the legend. Hovering over a detection displays its key parameters (e.g. its position and magnitude), and clicking a detection will search for the object in SIMBAD.
# 
# Detections are corrected for proper motion where possible, as stated by the ``corrected`` flag. Non-gaia detections are also corrected by "piggybacking" them with corresponding Gaia detections. The GALEX data that was found for van Maanen's Star from the :doc:`previous section <../getting_started/data_query>` is a match!
# 
# |
# |
#
# Overlay Configuration
# =====================
# By default, a detection overlay uses whichever band is listed first in the catalogue's overlay definition (see :doc:`here <../extension/overlays>`). If a different band is required, or perhaps multiple bands from the same catalogue need to be overlayed simultaneously, ``overlays`` can instead be passed as a ``dict``:

ps_query = query(
    "image",
    targets=2552928187080872832,
    survey="panstarrs",
    band="g",
    size=120,
    overlays={"galex": ["NUVmag", "FUVmag"]},
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
#    Where possible, detection markers are scaled in size by magnitude. Markers labelled ``galex detection`` in the above image are detections that had an invalid magnitude, and hence could not be scaled.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
