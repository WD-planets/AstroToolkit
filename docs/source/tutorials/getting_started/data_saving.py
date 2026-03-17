"""
##################
Local Data Storage
##################

To avoid having to rerun a query every time a script is used, :class:`DataSets <ATK.Models.DataSet>` of any kind can be stored as local fits files. These files can be read by ATK to re-create the original :class:`~ATK.Models.DataSet`.

|
|

Storing Data Automatically
==========================
:func:`~ATK.Tools.query` can be provided with a ``path``:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
# sphinx_gallery_end_ignore

from ATK import query

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex", path="example_data.fits")
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# This will automatically save the returned :class:`~ATK.Models.DataSet` to ``path``. Running the script again will read this local file instead of rerunning :func:`~ATK.Tools.query`.
#
# .. note::
#
#    Changes to any query parameters will automatically trigger the query to run again, overwriting the local file with the updated :class:`~ATK.Models.DataSet`.
#
# |
# |
#
# Storing Data Manually
# =====================
# A :class:`~ATK.Models.DataSet` can also be stored manually:

from ATK import query

target = 2552928187080872832
galex_query = query("vizier", targets=target, survey="galex")
galex_query.store("example_data.fits")
# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# If stored this way, the query will rerun every time the script is executed. The original :class:`~ATK.Models.DataSet` can be re-created with the :func:`~ATK.Tools.read` tool:

from ATK import read

data = read("example_data.fits")
data.show()
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
