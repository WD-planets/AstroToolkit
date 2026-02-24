"""
####################
Multi-Target Queries
####################

Queries are not limited to a single target. To perform a multi-target :func:`~ATK.Tools.query`, all we need to do is pass multiple valid targets (see the :doc:`first tutorial <data_query>` for examples of how to target a star).

Below is an example of a multi-target Vizier query for van Maanen's star and the cataclysmic variable Hu Leo:
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
    ctnr.data = df[df.columns[:7]]
# sphinx_gallery_end_ignore
galex_query.show()

# %%
# 
# | 
# |
#
# Accessing the Returned Data
# ===========================

# %%
# In the returned :class:`~ATK.Models.DataSet`, we now have two :class:`~ATK.Models.VizierEntry` containers. We could extract these as we did before:

van_maanen, hu_leo = galex_query.data

# %%
# 
# However, this can easily produce unexpected results if one or more of our queries do not return any data. Instead, we should use one of the fetch methods of the DataSet.
# 
# |
#
# Matching Returned Data by ID
# ----------------------------
# Since we used Gaia source IDs to target our stars, we can use the :meth:`~ATK.Models.DataSet.fetch_by_id` method, which returns a list of containers with a matching source ID:

van_maanen = galex_query.fetch_by_id(2552928187080872832)
hu_leo = galex_query.fetch_by_id(587316166180416640)

# sphinx_gallery_start_ignore
print(van_maanen)
print(hu_leo)
# sphinx_gallery_end_ignore

# %%
# 
# |
#
# Matching Returned Data by Coordinates
# -------------------------------------
# The equivalent method for searching by coordinates, :meth:`~ATK.Models.DataSet.fetch_by_coord`, which returns a list of containers with input positions (i.e. those that were provided to :func:`~ATK.Tools.query`) that fall within a given ``radius``:

from astropy.coordinates import SkyCoord
import astropy.units as u

coord = SkyCoord(12.2967, 5.3766, unit="deg", frame="icrs")
van_maanen = galex_query.fetch_by_coord(coord, radius=3*u.arcsec)

# sphinx_gallery_start_ignore
van_maanen
# sphinx_gallery_end_ignore

# %%
# 
# |
#
# Matching Returned Data by Target 
# --------------------------------
# Finally, :meth:`~ATK.Models.DataSet.fetch_by_target` can be used to retrieve containers using a :class:`~ATK.Models.Target` directly - allowing for exact matching without needing a ``radius``:

from ATK.Models import Target

coord_1 = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")
coord_2 = SkyCoord(141.1853, 8.0308, unit="deg", frame="icrs")

targets = [Target.from_coord(coord_1), Target.from_coord(coord_2)]

galex_query = query("vizier", targets=targets, survey="galex")

van_maanen = galex_query.fetch_by_target(targets[0])
hu_leo = galex_query.fetch_by_target(targets[1])

# sphinx_gallery_start_ignore
print(van_maanen)
print(hu_leo)
# sphinx_gallery_end_ignore

# %%
# .. note ::
#
#    As we have not used a Gaia Source ID to target van Maanen's Star here, no proper motion correction has been performed and hence no data was returned.
#
# |
# |
# |
#
# Download this Tutorial
# ======================
