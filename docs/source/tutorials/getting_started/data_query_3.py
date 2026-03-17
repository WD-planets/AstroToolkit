"""
####################
Multi-Target Queries
####################
Below is an example of a multi-target Vizier :func:`~ATK.Tools.query` for two sources: **van Maanen's star**, a lone white dwarf, and **Hu Leo**, a cataclysmic variable:
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
galex_query.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# 
# | 
# |
#
# Accessing the Returned Data
# ===========================
# The returned :class:`~ATK.Models.DataSet` contains two :class:`records <ATK.Models.Record>`. We could extract data from these as in the :doc:`previous tutorial <data_query>`:

van_maanen, hu_leo = galex_query.data

# %%
# 
# But this can easily produce unexpected results if one or more of our queries do not return any data. Instead, we should use one of the fetch methods of the :class:`~ATK.Models.DataSet`.
# 
# |
#
# Extracting Data by ID
# ---------------------
# Since Gaia source IDs were used in the query, the :meth:`~ATK.Models.DataSet.fetch_by_id` method can be used. This returns a list of containers with a matching source ID:

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
# Extracting Data by Coordinates
# ------------------------------
# The equivalent method for searching by coordinates is :meth:`~ATK.Models.DataSet.fetch_by_coord`, which returns a list of containers with **input** positions  - i.e. those that were provided to :func:`~ATK.Tools.query` - that fall within a given ``radius``:

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
# Finally, :meth:`~ATK.Models.DataSet.fetch_by_target` can be used to retrieve containers using a :class:`~ATK.Models.Target` directly. This allows for exact matching without using a ``radius``:

from ATK.Models import Target

coord_1 = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs") # van Maanen's Star
coord_2 = SkyCoord(141.1853, 8.0308, unit="deg", frame="icrs") # Hu Leo

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
#    As a Gaia Source ID was not used to target van Maanen's Star, no proper motion correction was performed and no data was returned.
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
