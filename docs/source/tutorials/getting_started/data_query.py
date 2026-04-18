"""
.. role:: bolditalic

   :class: bolditalic

#########################
Performing a Vizier Query
#########################
The most simple example of a query in ATK is a `Vizier <https://vizier.cds.unistra.fr/>`_ query. In this tutorial, data will be retrieved from GALEX. Queries of any kind use the :func:`~ATK.Tools.query` tool, which can be imported as follows:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from ATK.Config import SURVEY_CONFIG
from ATK.queries.vizier._query_info import SUPPORTED_SURVEYS

SURVEY_CONFIG.reset()

with open("../../auto_tutorials/getting_started/supported_aliases.rst", "w") as f:
    for alias in SUPPORTED_SURVEYS:
        f.write(f"    - {alias}\n")
# sphinx_gallery_end_ignore

from ATK import query

# %%
# or:

from ATK.Tools import query

# %%
# 
# |
# |
#
# Setting a Target
# ================
# The first step in performing a query is to set a target, which can be done in a few different ways. In this example we will use `van Maanen's Star <https://simbad.u-strasbg.fr/simbad/sim-id?Ident=van+Maanen%27s+Star&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id>`_.
#
# 1. Stars with Gaia DR3 data can be targeted by their Source ID:

target = 2552928187080872832

# %%
#
# |
#
# 2. The target's coordinates can be used directly by creating an astropy :class:`~astropy.coordinates.SkyCoord`:

from astropy.coordinates import SkyCoord

target = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")

# %%
# 
# |
# 
# 3. A :class:`~ATK.Models.Target` object can be created explicitly, either from a Gaia source ID or from the star's coordinates:

from ATK.Models import Target

target = Target.from_id(2552928187080872832)

# %% or:

from ATK.Models import Target
from astropy.coordinates import SkyCoord

coord = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")
target = Target.from_coord(coord)

# %%
# .. note::
#
#    Passing a :class:`~ATK.Models.Target` explicitly may seem unnecessary at first, but it offers a key advantage in queries that simultaneously target :doc:`multiple stars <../multiple_targets/multi_target_query>`.

# %%
# 
# | 
# |
# 
# Performing a Query
# ====================
# With a target set (**in this case via the star's coordinates**), we can now perform a GALEX `Vizier <https://vizier.cds.unistra.fr/>`_ query:

target = SkyCoord(12.2912, 5.3886, unit="deg", frame="icrs")
galex_query = query("vizier", targets=target, survey="galex")

# %%
# .. note::
#
#    By default, ATK supports queries via aliases to the following `Vizier <https://vizier.cds.unistra.fr/>`_ catalogues:
#
#    .. include:: supported_aliases.rst
# 
#    To utilise any other :doc:`Vizier <../extension/vizier>` catalogue, enter its catalogue ID (e.g. ``"I/355/gaiadr3"`` for Gaia DR3). A full tutorial on extending :func:`~ATK.Tools.query` to work with any `Vizier <https://vizier.cds.unistra.fr/>`_ catalogue (including automatic proper motion correction) can be found :doc:`here <../extension/vizier>`.
# 
# |
#
# The :func:`~ATK.Tools.query` tool returns a :class:`~ATK.Models.DataSet`, which contains key details of the request along with any returned data. To see the structure of the returned :class:`DataSet <ATK.Models.DataSet>`, we can call its :meth:`~ATK.Models.DataSet.show` method, which prints any ATK object in a human-readable format.
#
# Here, ``show_types=True`` is passed to :meth:`~ATK.Models.DataSet.show`. This forces the printing of all attribute data types:

galex_query.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# 
# The majority of the returned :class:`~ATK.Models.DataSet`'s attributes are immediately clear, but there are a few specific details to note:
# 
# - :attr:`~ATK.Models.DataSet.targets` lists the :class:`Targets <ATK.Models.Target>` of the search **as they were entered** (coordinates, frame and epoch).
#
# - :attr:`~ATK.Models.DataSet.exception` states whether any exceptions were encountered. If :attr:`~ATK.Models.DataSet.exception` is ``True``, then the query encountered something unexpected. The most common cause of this is that service being utilised (in this case `Vizier <https://vizier.cds.unistra.fr/>`_) is experiencing downtime.
#
# |
# 
# Any returned data is stored in the :class:`~ATK.Models.DataSet`'s :attr:`~ATK.Models.DataSet.data` attribute. Unfortunately, targeting van Maanen's Star with its coordinates has returned an empty :class:`~ATK.Models.DataSet`. We can instead try targeting it via its Gaia Source ID:

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex")
galex_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# **Targeting van Maanen's Star in this way has allowed ATK to utilise Gaia's astrometry to correct the position of our search for proper motion** :bolditalic:`before` **it is executed. Without increasing the query radius, we have retrieved GALEX data for van Maanen's Star.**
# 
# |
# |
#
# Using the Returned Data
# =======================
# All :class:`DataSets <ATK.Models.DataSet>` store any returned data as a list of containers matching the ``kind`` of data that was requested. Since we performed a `Vizier <https://vizier.cds.unistra.fr/>`_ query for one target, we get a list containing a single :class:`~ATK.Models.Record`:

galex_entry = galex_query.data[0]
# sphinx_gallery_start_ignore
print(galex_entry)
# sphinx_gallery_end_ignore

# %%
# As with any ATK object, :meth:`~ATK.Models.Record.show` can be used to see its structure:

galex_entry.show(show_types=True)
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# While the attributes of a :class:`~ATK.Models.DataSet` describe the query **as it was requested**, the attributes of any stored containers describe the query **as it was executed**:
#
# - :attr:`~ATK.Models.Record.catalogue` shows the ID of the requested :attr:`~ATK.Models.Record.survey` in Vizier.
#
# - :attr:`~ATK.Models.Record.correction` states the degree of proper motion correction that was achieved. A ``full`` correction indicates that the system has valid Gaia proper motion and distance, and so a complete 3-dimensional correction was performed. A ``partial`` correction occurs when the star has an invalid distance and so correction is purely angular on the sky - this is still fine in essentially all cases. If :attr:`~ATK.Models.Record.correction` is ``none``, no correction was performed - either because the star has invalid proper motion in Gaia, or because ATK does not know the median epoch of the requested survey (more on this later).
#
# - :attr:`~ATK.Models.Record.search_pos` gives the actual position (coordinates, frame and epoch) of the search. In this case, the coordinates of the target have been corrected from Gaia's epoch of January 2016 to GALEX's median epoch of August 2006.
# 
# |
#
# A :class:`~ATK.Models.Record` stores its data as an astropy :class:`~astropy.table.Table`. Any parameter can therefore be extracted from a catalogue (in this case the GALEX Near-UV magnitude) as follows:

nuv_mag = galex_entry.table["NUVmag"][0]

# sphinx_gallery_start_ignore
nuv_mag
# sphinx_gallery_end_ignore

# %%
# 
# |
# 
# .. note::
#    
#    `Vizier <https://vizier.cds.unistra.fr/>`_ queries to Gaia will always return exactly one row when targeting a star via a valid Gaia Source ID. Any other search combination will return all rows within the query radius.
#
# |
# |
#
# Configuring a Query 
# ===================
# The arguments that can be passed to :func:`~ATK.Tools.query` depend on the ``kind`` of data that is being requested. `Vizier <https://vizier.cds.unistra.fr/>`_ queries do not have many optional arguments, but one example is the ``radius``, which defaulted to 3 arcseconds above. The best way to override this is to set the ``radius`` explicity using :mod:`astropy.units`:

from ATK import query
import astropy.units as u

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex", radius=2*u.arcmin)
galex_query.show()
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
#
# **This automatically converts all structure attributes to match the requested spatial scale (in this case,** :attr:`~ATK.Models.Record.separation` **is now given in arcminutes).**
#
# |
#
# .. note::
#
#    While a default radius of 3 arcseconds is reasonable for a modern optical telescope, for other workloads it may be preferrable to work at a different spatial resolution. To accomodate this, most of the core defaults in ATK are configurable via a config file. See :doc:`here <../configuration/config>` for details.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
