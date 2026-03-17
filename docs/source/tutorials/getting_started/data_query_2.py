"""
###############################
Vizier Queries to Any Catalogue
###############################
ATK uses aliases to map between Vizier catalogue names (``"galex"``) and Vizier catalogue IDs (``"II/335/galex_ais"``). In this tutorial, :func:`~ATK.Tools.query` will be extended to work with AllWISE.

|

To immediately access any Vizier catalogue, its catalogue ID can be passed directly to :func:`~ATK.Tools.query`:
"""

# sphinx_gallery_start_ignore
from ATK.Config import ALIAS_CONFIG, EPOCH_CONFIG
from ATK.queries.vizier._query_info import SUPPORTED_SURVEYS

ALIAS_CONFIG.reset()
EPOCH_CONFIG.reset()

with open("../../auto_tutorials/getting_started/supported_aliases.rst", "w") as f:
    for alias in SUPPORTED_SURVEYS:
        f.write(f"    - {alias}\n")
# sphinx_gallery_end_ignore

from ATK import query

target = 2552928187080872832
allwise_query = query("vizier", targets=target, survey="II/328/allwise")
allwise_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
#
# .. note::
#
#    By default, ATK supports the following aliases to Vizier catalogues:
#
#    .. include:: supported_aliases.rst
#
# |
#
# **However, if a Vizier catalogue is being used regularly it may become tedious to remember its ID. Additionally, and perhaps more importantly, we are not benefitting from proper motion correction as ATK doesn't know AllWISE's epoch.**
#
# |
# |
#
# Adding a New Catalogue Alias
# ============================
# To solve the first of these problems, we can add a new catalogue alias. This is done through an alias config file, which can be accessed from the command line:
#
# .. code-block:: console
#
#    $ ATKalias show
#
# or from inside a script:

from ATK.Config import ALIAS_CONFIG

ALIAS_CONFIG.show()

# %%
# |
#
# Editing the Alias File
# ----------------------
# The easiest way to edit the alias file is to open it in the default text editor:
#
# .. code-block:: console
#
#    $ ATKalias open
#
# .. code-block:: python
#
#    from ATK.Config import ALIAS_CONFIG
#
#    ALIAS_CONFIG.open()
#
# |
#
# Utilising an Alias
# ------------------
# After adding ``allwise = II/328/allwise`` to the alias file, ``allwise`` can now be used as an alias in Vizier queries:

# sphinx_gallery_start_ignore
ALIAS_CONFIG.vizier_aliases.allwise = "II/328/allwise"
# sphinx_gallery_end_ignore

allwise_query = query("vizier", targets=target, survey="allwise")
allwise_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore


# %%
#
# |
# |
#
# Setting Survey Epochs
# =====================
# An alias to AllWISE has now been added, but we are still getting no data. This is because ATK does not know AllWISE's epoch, and so no correction can be performed. ATK stores the epochs of all surveys in an epoch file, which can be accessed in the same way as the alias file:
#
# .. code-block:: console
#
#    $ ATKepoch show

from ATK.Config import EPOCH_CONFIG

EPOCH_CONFIG.show()

# %%
# |
#
# Editing the Epoch File
# ----------------------
# Just like the alias file, the epoch file can be opened in the default text editor with:
#
# .. code-block:: console
#
#    $ ATKepoch open
#
# .. code-block:: python
#
#    from ATK.Config import EPOCH_CONFIG
#
#    EPOCH_CONFIG.open()
#
# An epoch can then be set for AllWISE by adding ``allwise = 2010-06-01T00:00:00.000`` under the category ``vizier_aliases``.
#
# |
#
# Rerunning the Query
# -------------------
# Running the same query as above now utilises proper motion correction to retrieve AllWISE data for van Maanen's Star:

# sphinx_gallery_start_ignore
EPOCH_CONFIG.vizier_aliases.allwise = "2010-05-01T00:00:00.000"
# sphinx_gallery_end_ignore

target = 2552928187080872832
allwise_query = query("vizier", targets=target, survey="allwise")
allwise_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
