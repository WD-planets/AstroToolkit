"""
##############
Exporting Data
##############
If data is needed for external use, it may be helpful to export ATK data into a tabular format (i.e. a :class:`~pandas.DataFrame` or :class:`~astropy.table.Table`).

This can be done via the **IO methods** :meth:`~ATK.base.DataFrameIOMixin.to_dataframe` and :meth:`~ATK.base.TableIOMixin.to_table`.

|

.. note::

   The availability of IO methods depends on the kind of data that is being used. See the documentation of the specific data container for information (e.g. :class:`here <ATK.Models.Spectrum>` if working with :class:`~ATK.Models.Spectrum` objects).

|

In this tutorial, the aforementioned **IO methods** will be demonstrated on a :class:`~ATK.Models.Spectrum`:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from IPython.display import display
# sphinx_gallery_end_ignore

from ATK import query

sdss_query = query("spectrum", targets=2802295461460467072, survey="sdss")
sdss_query.show()

# %%
# |
#  
# To convert the stored :class:`~ATK.Models.Spectrum` into a :class:`~pandas.DataFrame`, simply extract it from the :class:`~ATK.Models.DataSet` and call :meth:`~ATK.base.DataFrameIOMixin.to_dataframe`:

df = sdss_query.data[0].to_dataframe()
# sphinx_gallery_start_ignore
df
# sphinx_gallery_end_ignore

# %%
# or to convert to a :class:`~astropy.table.Table`:

table = sdss_query.data[0].to_table()
# sphinx_gallery_start_ignore
table
# sphinx_gallery_end_ignore

# %%
# |
#
# .. note::
#
#    In either case, non-array-like attributes are not included in the output :class:`~pandas.DataFrame` or :class:`~astropy.table.Table`. If required, these should be extracted from the Spectrum alongside the use of :meth:`~ATK.base.DataFrameIOMixin.to_dataframe` or :meth:`~ATK.base.TableIOMixin.to_table`.
#
#    Also note that :meth:`~ATK.base.TableIOMixin.to_table` offers a key advantage over :meth:`~ATK.base.DataFrameIOMixin.to_dataframe`, in that units are retained and stored in the output :class:`~astropy.table.Table`.

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
