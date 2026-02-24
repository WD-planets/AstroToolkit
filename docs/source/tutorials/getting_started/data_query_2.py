"""
##########################
Generalised Vizier Queries
##########################
To immediately use any Vizier catalogue, we can simply pass the ``catalogue`` parameter with a Vizier catalogue ID instead of specifying a ``survey``. For this tutorial, we will extend our query to also work for AllWISE:
"""

from ATK import query

target = 2552928187080872832
allwise_query = query("vizier", targets=target, catalogue="II/328/allwise")
allwise_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# This works, but if we regularly use AllWISE it may become tedious to remember its Vizier ID. Additionally (and perhaps more importantly), we are not benefitting from proper motion correction as ATK doesn't know AllWISE's median epoch.
#
# |
# |
#
# Adding a New Vizier Catalogue
# =============================
# To solve the first of these problems, we can add a new catalogue alias. ATK uses an alias file to map between surveys (AllWISE) and the corresponding Vizier catalogue ID (II/328/allwise). We can see the current state of the alias file from the command line:
#
# .. code-block:: console
#
#    $ ATKalias show
#
# or from inside a script:

from ATK.Config import ALIAS_CONFIG

ALIAS_CONFIG.show()

# %%
# The easiest way to edit the alias file is to simply open it in the default text editor with:
#
# .. code-block:: console
#
#    $ ATKalias open
#
# from the command line, or:
#
# .. code-block:: python
#
#    ALIAS_CONFIG.open()
#
# from inside a script.
#
# |
#
# After adding the following line to the alias file:
#
# .. code-block:: console
#
#    allwise = II/328/allwise
#
# we can now use ``allwise`` as an alias inside Vizier queries:

# sphinx_gallery_start_ignore
ALIAS_CONFIG.vizier_aliases.allwise = "II/328/allwise"
# sphinx_gallery_end_ignore

allwise_query = query("vizier", targets=target, survey="allwise")
allwise_query.show()


# %%
#
# |
# |
#
# Setting a Catalogue's Epoch
# ===========================
# We have now added an catalogue alias to AllWISE, but we are still getting no data as we have not yet set AllWISE's epoch. ATK stores the epochs of all supported surveys in an epoch file, which (like the alias file) can be accessed via the command line or from inside a script:
#
# .. code-block:: console
#
#    $ ATKepoch show
#
# or:

from ATK.Config import EPOCH_CONFIG

EPOCH_CONFIG.show()

# %%
#
# We can then use:
#
# .. code-block:: console
#
#    $ ATKepoch open
#
# or:
#
# .. code-block:: python
#
#    from ATK.Config import EPOCH_CONFIG
#
#    EPOCH_CONFIG.open()
#
# to add the median epoch of AllWISE (mid-2010) under the category **vizier_aliases**:
#
# .. code-block:: console
#
#    allwise = 2010-06-01T00:00:00.000
#
# |
# |
#
# If we now rerun our query, we successfully retrieve AllWISE data for van Maanen's Star:

# sphinx_gallery_start_ignore
EPOCH_CONFIG.vizier_aliases.allwise = "2010-05-01T00:00:00.000"
# sphinx_gallery_end_ignore

target = 2552928187080872832
allwise_query = query("vizier", targets=target, catalogue="allwise")
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
# Download this Tutorial
# ======================
