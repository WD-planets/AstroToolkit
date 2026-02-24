"""
#########################
Performing a Vizier Query
#########################
Before moving onto plottable queries, we will start with the most simple example: a `Vizier <https://vizier.cds.unistra.fr/>`_ query - in this case to retrieve data from GALEX. All queries use the :func:`~ATK.Tools.query` tool, which can be imported as follows:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
# sphinx_gallery_end_ignore

from ATK import query

# %%
# 
# |
# |
#
# Setting a Target
# ================
#
# The first step in performing a query is to set our target, which can be done in a few different ways. In this example we will again use `van Maanen's Star <https://simbad.u-strasbg.fr/simbad/sim-id?Ident=van+Maanen%27s+Star&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id>`_.
#
# 1. We can target our star by using its Gaia Source ID:

target = 2552928187080872832

# %%
# 2. We can provide a :class:`~astropy.coordinates.SkyCoord` with the target's coordinates:

from astropy.coordinates import SkyCoord

target = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")

# %%
# 3. We can explicitly create a :class:`~ATK.Models.Target` object, either from a Gaia source ID or from the star's coordinates in the form of a SkyCoord:

from ATK.Models import Target

target = Target.from_id(2552928187080872832)

# %% or:

from ATK.Models import Target
from astropy.coordinates import SkyCoord

coord = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")
target = Target.from_coord(coord)

# %%
# 
# | 
# |
# 
# Performing the Query
# ====================
# With our target set (in this case via the star's coordinates), we can now fetch some data from Vizier:

target = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")
galex_query = query("vizier", targets=target, survey="galex")

# sphinx_gallery_start_ignore
from ATK.queries.vizier._query_info import SUPPORTED_SURVEYS
with open("../../auto_tutorials/getting_started/supported_aliases.rst", "w") as f:
    f.write(".. note::\n")
    f.write("    By default, ATK supports the following aliases to Vizier catalogues:\n")
    for alias in SUPPORTED_SURVEYS:
        f.write(f"        - {alias}\n")
    f.write("\n")
    f.write("    Instructions on utilising any other Vizier catalogues can be found in the :doc:`next tutorial <data_query_2>`.")
# sphinx_gallery_end_ignore

# %%
# 
# .. include:: supported_aliases.rst 
#
# The :func:`~ATK.Tools.query` tool returns a :class:`~ATK.Models.DataSet` which contains information about our query and any returned data. We can see the structure of the returned :class:`DataSet <ATK.Models.DataSet>` by calling its :meth:`~ATK.Models.DataSet.show` method:

galex_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# We can see at the top that our query returned a :class:`~ATK.Models.DataSet` of GALEX Vizier data, and the attributes of this :class:`~ATK.Models.DataSet` are shown underneath (the first two of which are self-explanatory). Following this, the :class:`~ATK.Models.DataSet` has three additional attributes that tell us information about the query we have requested:
#
# - The :attr:`~ATK.Models.DataSet.targets` attribute shows the targets that we entered for the search (coordinates, frame and epoch).
#
# - The :attr:`~ATK.Models.DataSet.radius` attribute gives the radius of the cone (in this case three arcseconds).
#
# - The :attr:`~ATK.Models.DataSet.exception` attribute tells us that no exceptions were encountered. If :attr:`~ATK.Models.DataSet.exception` were instead ``True``, this would suggest that the API that we are utilising (in this case Vizier) is experiencing issues.
#
# Finally, the :attr:`~ATK.Models.DataSet.data` attribute stores any returned data. Unfortunately, targeting the system with its coordinates (taken from `SIMBAD <https://simbad.cds.unistra.fr/simbad/>`_) has not returned any data. We can instead try to target it with its Gaia source ID (as above):

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex")
galex_query.show()

# %%
# The key difference here is that targeting our star via its Gaia source ID has allowed us to use Gaia's astrometry to correct the position of our search before it is executed. With the same query radius, we have now retrieved GALEX data for van Maanen's Star!
# 
# | 
# |
#
# Using the Returned Data
# =======================
# Now that we have successfully acquired GALEX data, we can extract it from our :class:`~ATK.Models.DataSet`. All returned data in ATK is stored as a list of containers, with the exact container type depending on the kind of data that we requested. In this case we get a list containing a single :class:`~ATK.Models.VizierEntry`:

galex_entry = galex_query.data[0]

# %%
# As with any ATK data structure, we can use the :meth:`~ATK.Models.VizierEntry.show` method to see the structure of the :class:`~ATK.Models.VizierEntry` (we could actually already see it nested above when we called :meth:`~ATK.Models.DataSet.show` on the returned :class:`~ATK.Models.DataSet`):

galex_entry.show()

# %%
# Just as the attributes of the returned :class:`~ATK.Models.DataSet` described the query that we requested, those of the :class:`~ATK.Models.VizierEntry` provide the details of the final query as it was performed:
#
# - The :attr:`~ATK.Models.VizierEntry.survey` attribute gives the name of the survey to which the data pertains.
#
# - The :attr:`~ATK.Models.VizierEntry.catalogue` attribute gives the ID of this survey in Vizier.
#
# - The :attr:`~ATK.Models.VizierEntry.correction` attribute tells us the degree of proper motion correction that was achieved. A ``full`` correction indicates that the system has valid Gaia proper motion and distance, and so a complete 3-dimensional correction was performed. A ``partial`` correction occurs when the star has an invalid distance and so correction is purely angular on the sky - this is still fine in essentially all cases. A :attr:`~ATK.Models.VizierEntry.correction` of ``none`` means that no correction was performed, either because the system has invalid proper motion in Gaia, or because an epoch has not been set for the requested survey (more on this later).
#
# - The :attr:`~ATK.Models.VizierEntry.search_pos` attribute gives us the actual position (coordinates, frame and epoch) on the sky at which the cone search was performed (notice how the coordinates and epoch differ from the attributes of the data set due to the correction that was performed).
#
# - The :attr:`~ATK.Models.VizierEntry.separation` attribute gives the separation on the sky between the position of the cone search and the returned data.

# %%
# We can now get the returned :class:`~pandas.DataFrame` and hence utilise any parameter from the Vizier catalogue, e.g. the Near-UV (NUV) magnitude:

nuv_mag = galex_entry.data["NUVmag"][0]

# sphinx_gallery_start_ignore
nuv_mag
# sphinx_gallery_end_ignore

# %%
# 
# |
# |
#
# Further Configuration
# =====================
# The arguments that can be passed to :func:`~ATK.Tools.query` depend on the ``kind`` of data that is being requested. Vizier queries do not have many optional arguments, but one example is the ``radius`` of our query, which defaulted to 3 arcseconds above. The easiest way to override this is to set the ``radius`` explicity using :mod:`astropy.units`:

from ATK import query
import astropy.units as u

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex", radius=5*u.arcsec)
galex_query.show()

# %%
# .. note::
#
#    While a default radius of 3 arcseconds is reasonable for a modern optical telescope, for other workloads it may be preferrable to work at a different spatial resolution. To accomodate this, most of the core defaults in ATK are configurable via a config file. See :doc:`here <../configuration/config>` for details.

# %%
# |
# |
# |
#
# Download this Tutorial
# ======================
