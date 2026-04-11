"""
####################
Multi-Target Queries
####################
All tasks in ATK (including all of those showcased in previous tutorials) are designed to work modularly with multiple targets. **A target matching system is used to internally link each target (whether that be a Gaia source ID, a** :class:`~ATK.Models.Target` **or a** :class:`~astropy.coordinates.SkyCoord`) **to the data that it returns.**

|

Performing a Multi-Target Query
===============================
A multi-target query can be performed by passing multiple targets to :func:`~ATK.Tools.query`. An example of a multi-target `Vizier <https://vizier.cds.unistra.fr/>`_ query for two targets - **van Maanen's star** and **Hu Leo**, a cataclysmic variable - is shown below:
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
    tbl = ctnr.table
    ctnr.table = tbl[tbl.colnames[:5]]
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
# Accessing Multi-Target Data
# ===========================
# The returned :class:`~ATK.Models.DataSet` contains two :class:`Records <ATK.Models.Record>`. These could be extracted by pulling them out of the :attr:`~ATK.Models.DataSet.data` attribute:

van_maanen, hu_leo = galex_query.data

# %%
#
# **But this can easily produce unexpected results if no data is returned for one or more of the targets.** Instead, the :class:`~ATK.Models.DataSet`'s :meth:`~ATK.Models.DataSet.split` method should be used. :meth:`~ATK.Models.DataSet.split` takes any valid target (i.e. a Gaia source ID, a :class:`~ATK.Models.Target`, or an astropy :class:`~astropy.coordinates.SkyCoord`) and returns a split of the :class:`~ATK.Models.DataSet` which only contains that target's data.
#
# |
# 
# Splitting by ID 
# ---------------
# Since Gaia source IDs were used in this example, these can be passed to :meth:`~ATK.Models.DataSet.split`:

van_maanen_2 = galex_query.split(2552928187080872832, inplace=False)
van_maanen_2.show()

# %%
# **The returned** :class:`~ATK.Models.DataSet` **contains only the requested target in its** :attr:`~ATK.Models.DataSet.targets` **attribute, and only the corresponding** :class:`~ATK.Models.Record` **in its** :attr:`~ATK.Models.DataSet.data` **attribute.**
# 
# **Multiple targets can also be passed to** :meth:`~ATK.Models.DataSet.split` **(in this case, passing both targets to** :meth:`~ATK.Models.DataSet.split` **would leave the** :class:`~ATK.Models.DataSet` **unchanged.)**
# 
# .. note::
#
#    ``inplace = False`` tells :meth:`~ATK.Models.DataSet.split` to act on and return a copy of the :class:`~ATK.Models.DataSet` (leaving the original unchanged).
#
# |
#
# Splitting by Coordinates 
# ------------------------
# Equivalently, an astropy :class:`~astropy.coordinates.SkyCoord` can be passed to :meth:`~ATK.Models.DataSet.split`. This returns a split of the original :class:`~ATK.Models.DataSet` which only contains data corresponding to **input positions** (i.e. the uncorrected coordinates that were provided to :func:`~ATK.Tools.query` by the user) that fall within a given ``radius``:

import astropy.units as u
from astropy.coordinates import SkyCoord

coord = SkyCoord(141.185, 8.031, unit="deg", frame="icrs")
hu_leo = galex_query.split(coord, inplace=False, radius = 5 * u.arcsec)
hu_leo.show()

# %%
#
# .. note::
#
#    In queries that target stars with their Gaia source IDs (like the one performed here), splitting via a :class:`~astropy.coordinates.SkyCoord` isn't particularly useful. However, it is worth noting that the **input positions** of the stars in this case would be those given by Gaia DR3.
#
# |
#
# Splitting by Target 
# -------------------
# Finally, a :class:`~ATK.Models.Target` can be passed to :meth:`~ATK.Models.DataSet.split`. **This enables exact matching without the need of a** ``radius``, **and is therefore preferable if a Gaia source ID isn't available:**

from ATK.Models import Target

coord_1 = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")  # van Maanen's Star
coord_2 = SkyCoord(141.1853, 8.0308, unit="deg", frame="icrs")  # Hu Leo
targets = [Target.from_coord(coord_1), Target.from_coord(coord_2)]

gaia_query = query("vizier", targets=targets, survey="gaia")
# sphinx_gallery_start_ignore
for ctnr in gaia_query.data:
    tbl = ctnr.table
    ctnr.table = tbl[tbl.colnames[:9]]
# sphinx_gallery_end_ignore

hu_leo = gaia_query.split(targets[1], inplace=False)
hu_leo.show()

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
