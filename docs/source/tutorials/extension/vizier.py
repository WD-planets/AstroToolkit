"""
#########################
Setting Survey Parameters
#########################
The Survey File
===============
ATK stores all survey parameters in a file, which can be opened either from the commandline:

.. code-block:: console

   $ ATKsurvey open

or from inside a script:

.. code-block:: python

   from ATK.Config import SURVEY_CONFIG

   SURVEY_CONFIG.open()

|

To see the current state of the **survey file**, use:

.. code-block:: console

   $ ATKsurvey show
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK.Config import SURVEY_CONFIG

# sphinx_gallery_start_ignore
SURVEY_CONFIG.reset()
# sphinx_gallery_end_ignore
SURVEY_CONFIG.show()

# %%
# |
# |
#
# Adding a Catalogue Alias
# ========================
# Any `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue can be immediately accessed by passing its catalogue ID to :func:`~ATK.Tools.query`:

from ATK import query

allwise_query = query("vizier", targets=2552928187080872832, survey="II/328/allwise")
# sphinx_gallery_start_ignore
pass 
# sphinx_gallery_end_ignore


# %%
# |
# 
# However, it may become tedious to remember catalogue IDs in prolonged workloads. Additionally, and perhaps more importantly, **the above does not benefit from proper motion correction as ATK doesn't know the catalogue's epoch.** Both of these problems can be solved by adding a catalogue alias to the **survey file**.
#
# In this tutorial, :func:`~ATK.Tools.query` will be extended with **AllWISE**.
# 
# |
# 
# Editing the Survey File 
# -----------------------
# To add an **AllWISE** alias, the following must be added to the **survey file**, under ``vizier``:
#
# .. code-block:: console
#
#    allwise
#        id = II/328/allwise
# 
# |
# 
# This can be done manually by opening **survey file** in the default text editor, or by setting it directly:
#
# .. code-block:: console
# 
#    $ ATKsurvey set vizier allwise --id II/328/allwise
# 
# .. code-block:: python 
# 
#    from ATK.Config import SURVEY_CONFIG
# 
#    SURVEY_CONFIG["vizier"]["allwise"]["id"] = "II/328/allwise"
#
# |
# 
# Testing the New Alias 
# ---------------------

# sphinx_gallery_start_ignore
SURVEY_CONFIG["vizier"]["allwise"]["id"] = "II/328/allwise"
# sphinx_gallery_end_ignore
from ATK import query

target = 2552928187080872832
allwise_query = query("vizier", targets=target, survey="allwise", path="allwise_data_1.fits.gz")
allwise_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %% 
# |
# 
# The new alias is working, but **the query has not returned any data.** This is because **no epoch information has been added.** This can be amended by setting the ``epoch`` key with the median epoch of **AllWISE** (approximately ``2010-05-01T00:00:00.000``):
#
# .. code-block:: console
# 
#    $ ATKsurvey set vizier allwise --epoch 2010-05-01T00:00:00.000
# 
# .. code-block:: python 
# 
#    from ATK.Config import SURVEY_CONFIG
# 
#    SURVEY_CONFIG["vizier"]["allwise"]["epoch"] = "2010-05-01T00:00:00.000"

# %% 
# | 
# 
# **Upon running the same query again, data is now returned:**

# sphinx_gallery_start_ignore
SURVEY_CONFIG["vizier"]["allwise"]["epoch"] = "2010-05-01T00:00:00.000"
# sphinx_gallery_end_ignore
from ATK import query

target = 2552928187080872832
allwise_query = query("vizier", targets=target, survey="allwise", path="allwise_data_2.fits.gz")
allwise_query.show()

# %%
#
# |
# |
#
# Defining an Image Overlay
# =========================
# The **AllWISE** `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue alias added above can be extended into an image overlay definition by adding some extra information. **As a minimum, all overlay definitions require three columns for positional information: a longitude column, a latitude column, and a frame on which these are based.**
#
# In the case of **AllWISE**, this means adding the following:
#
# .. code-block:: console
#
#    $ ATKsurvey set vizier allwise --lon RAJ2000 --lat DEJ2000 --frame icrs
#
# .. code-block:: python
# 
#    from ATK.Config import SURVEY_CONFIG
#
#    SURVEY_CONFIG["vizier"]["allwise"]["lon"] = "RAJ2000"
#    SURVEY_CONFIG["vizier"]["allwise"]["lat"] = "DEJ2000"
#    SURVEY_CONFIG["vizier"]["allwise"]["frame"] = "icrs"
# 
# | 
# 
# Once these keys have been added, an **AllWISE** overlay can be requested in an image :func:`~ATK.Tools.query`:

# sphinx_gallery_start_ignore
SURVEY_CONFIG["vizier"]["allwise"]["lon"] = "RAJ2000"
SURVEY_CONFIG["vizier"]["allwise"]["lat"] = "DEJ2000"
SURVEY_CONFIG["vizier"]["allwise"]["frame"] = "icrs"
# sphinx_gallery_end_ignore
from ATK import query
import astropy.units as u

target = 2552928187080872832
panstarrs_query = query("image", targets=target, survey="panstarrs", band="g", size=90*u.arcsec, overlays=["allwise"], path="image_data_1.fits.gz")
# sphinx_gallery_start_ignore
panstarrs_query.plot()
figure = format_plot(panstarrs_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
panstarrs_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Adding Photometric Information
# ------------------------------
# Providing positional information is already enough to generate corrected detections - but for photometric surveys like **AllWISE**, magnitude information can also be added by supplying a list of **magnitude** and **magnitude error** columns:
#
# .. code-block:: console
#
#    $ ATKsurvey set vizier allwise --mags W1mag W2mag W3mag W4mag --errors e_W1mag e_W2mag e_W3mag e_W4mag
#
# .. code-block:: python
# 
#    from ATK.Config import SURVEY_CONFIG
#
#    SURVEY_CONFIG["vizier"]["allwise"]["mags"] = ["W1mag", "W2mag", "W3mag", "W4mag"]
#    SURVEY_CONFIG["vizier"]["allwise"]["errors"] = ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"]
# 
# |
# 
# This turns the **AllWISE** detection overlay into a **photometric** detection overlay, which scales marker sizes with magnitude and adds photometric information to the hover interaction:

# sphinx_gallery_start_ignore
SURVEY_CONFIG["vizier"]["allwise"]["mags"] = ["W1mag", "W2mag", "W3mag", "W4mag"]
SURVEY_CONFIG["vizier"]["allwise"]["errors"] = ["e_W1mag", "e_W2mag", "e_W3mag", "e_W4mag"]
# sphinx_gallery_end_ignore
target = 2552928187080872832
panstarrs_query = query("image", targets=target, survey="panstarrs", band="g", size=90*u.arcsec, overlays=["allwise"], path="image_data_2.fits.gz")
# sphinx_gallery_start_ignore
panstarrs_query.plot()
figure = format_plot(panstarrs_query.figure, 1.5, 1.5, True)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
panstarrs_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
# |
#
# Additional Survey File Operations
# =================================
# Resetting the Survey File
# -------------------------
# To reset the **survey file** to its default state:
#
# .. code-block:: console
#
#    $ ATKsurvey reset
#
# .. code-block:: python
#
#    from ATK.Config import SURVEY_CONFIG
#
#    SURVEY_CONFIG.reset()
#
# |
#
# Deleting an Alias
# -----------------
# To delete an alias (in this case the **AllWISE** alias created in this tutorial) from the **survey file**, either delete it from the commandline:
#
# .. code-block:: console
#
#    $ ATKsurvey del allwise
#
# or set it to ``None``:
#
# .. code-block:: python
#
#    from ATK.Config import SURVEY_CONFIG
#
#    SURVEY_CONFIG["vizier"]["allwise"] = None

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
