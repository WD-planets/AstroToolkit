"""
##################
Detection Overlays
##################
To simplify the matching of positional data, :class:`Images <ATK.Models.Image>` support detection overlays.

.. note::

   By default, detection overlays are configured for the following Vizier catalogues:

   .. include:: ../getting_started/supported_aliases.rst

   A tutorial on defining detection overlays for any Vizier catalogue can be found :doc:`here <../extension/overlays>`.

|
|

Requesting an Overlay
=====================
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
