"""
####################
Multi-Target Queries
####################
All tasks in ATK are designed to work modularly with any number of targets. Below is an example of a multi-target `Vizier <https://vizier.cds.unistra.fr/>`_ query for two sources: **van Maanen's star** and **Hu Leo**, a cataclysmic variable:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
# sphinx_gallery_end_ignore

from ATK import query

targets = [2552928187080872832, 587316166180416640]
galex_query = query("vizier", targets=targets, survey="gaia")
# sphinx_gallery_start_ignore
for ctnr in galex_query.data:
    df = ctnr.data
    ctnr.data = df[df.columns[:5]]
# sphinx_gallery_end_ignore
galex_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# .. note::
# 
#    The returned :class:`Records <ATK.Models.Record>` have been truncated here for clarity.
# 
# | 
# |
#
# Accessing the Returned Data
# ===========================
# The returned :class:`~ATK.Models.DataSet` contains two :class:`Records <ATK.Models.Record>`. We could extract these as in the :doc:`previous tutorial <data_query>`:

van_maanen, hu_leo = galex_query.data

# %%
# 
# **But this can easily produce unexpected results if no data is returned for one or more of the targets.** Instead, one of the **fetch methods** of the :class:`~ATK.Models.DataSet` should be used.
# 
# |
#
# Extracting Data by ID
# ---------------------
# Since Gaia source IDs were used in the query, the :meth:`~ATK.Models.DataSet.fetch_by_id` method can be used. This returns a list of containers with a matching source ID:

van_maanen_2 = galex_query.fetch_by_id(2552928187080872832)[0]

van_maanen_2.show()

# %%
# 
# |
#
# Extracting Data by Coordinates
# ------------------------------
# The equivalent method for searching by coordinates is :meth:`~ATK.Models.DataSet.fetch_by_coord`, which returns a list of containers with **input positions**  - i.e. those that were provided to :func:`~ATK.Tools.query` - that fall within a given ``radius``:

from astropy.coordinates import SkyCoord
import astropy.units as u

coord = SkyCoord(141.185, 8.031, unit="deg", frame="icrs")
hu_leo = galex_query.fetch_by_coord(coord, radius=3*u.arcsec)[0]

hu_leo.show()

# %%
# 
# .. note::
# 
#    In queries that target stars with their Gaia source IDs (like the one performed here), using :meth:`~ATK.Models.DataSet.fetch_by_coord` likely doesn't make much sense. However, it is worth noting that the **input positions** of the stars in this case are those in Gaia DR3.
#
# |
#
# Matching Returned Data by Target 
# --------------------------------
# Finally, :meth:`~ATK.Models.DataSet.fetch_by_target` can be used to retrieve containers using a :class:`~ATK.Models.Target` directly. This allows for exact matching without the need of a ``radius``:

from ATK.Models import Target

coord_1 = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs") # van Maanen's Star
coord_2 = SkyCoord(141.1853, 8.0308, unit="deg", frame="icrs") # Hu Leo

targets = [Target.from_coord(coord_1), Target.from_coord(coord_2)]

gaia_query = query("vizier", targets=targets, survey="gaia")
# sphinx_gallery_start_ignore
for ctnr in gaia_query.data:
    df = ctnr.data
    ctnr.data = df[df.columns[:9]]
# sphinx_gallery_end_ignore

hu_leo = gaia_query.fetch_by_target(targets[1])[0]

hu_leo.show()

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
