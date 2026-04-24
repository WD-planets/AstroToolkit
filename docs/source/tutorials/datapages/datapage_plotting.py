"""
####################
Generating Datapages
####################
ATK supports the creation of **datapages** as a means of neatly combining plots into a single grid. To demonstrate this, some data must first be acquired:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from bokeh.document import Document
from ATK.Config import CONFIG
CONFIG.reset()
CONFIG.datapage_settings.grid_size = 150
CONFIG.datapage_settings.font_size = 6
# sphinx_gallery_end_ignore
from ATK import query

target = 587316166180416640

lc = query("lightcurve", survey="asassn", targets=target, path="datapage_lc.fits.gz")
pspec = lc.apply("pspec", fmin=0, fmax=60, samples=10000, inplace=False)
fold = lc.apply("fold", fmin=0, fmax=60, samples=10000, inplace=False)
image = query("image", survey="panstarrs", band="g", targets=target, overlays=["gaia","galex"], path="datapage_image.fits.gz")
spec = query("spectrum", survey="sdss", targets=target, path="datapage_spec.fits.gz")
sed = query("sed", targets=target, path="datapage_sed.fits.gz")
hrd = query("hrd", targets=target, path="datapage_hrd.fits.gz")
table = query("datatable", columns={"gaia": ["Gmag", "BPmag", "RPmag"], "galex": ["NUVmag", "FUVmag"]}, targets=target, radius=image.data[0].size)

# %%
# The only other requirement is to define the **datapage's** layout. The easiest way to do this is by passing a list of rows of :class:`DataSets <ATK.Models.DataSet>`:

layout = [[image, image, hrd,   hrd,   lc,    lc,    lc,    lc   ],
          [image, image, hrd,   hrd,   lc,    lc,    lc,    lc   ],
          [pspec, pspec, pspec, pspec, fold,  fold,  fold,  fold ],
          [pspec, pspec, pspec, pspec, fold,  fold,  fold,  fold ],
          [spec,  spec,  spec,  spec,  sed,   sed,   sed,   sed  ],
          [spec,  spec,  spec,  spec,  sed,   sed,   sed,   sed  ],
          [table, table, table, table, table, table, table, table],
          [table, table, table, table, table, table, table, table]]

# %%
# This defines a **datapage** where the first row contains a 2x2 image in the top-left corner, followed by a 2x2 HRD and a 4x2 light curve, etc. Passing this to :func:`~ATK.Visualisation.grid` returns a :class:`~ATK.Models.DataPages` object, which contains one datapage per target:

from ATK import grid

datapages = grid(layout)
datapages.show()

# %%
# |
# |
#
# Opening a Datapage 
# ==================
# The only remaining step is to open or save the **datapage**. Like everything in ATK, :class:`DataPages <ATK.Models.DataPages>` are designed to work with multiple targets. A **datapage** for a specific target can be opened by passing any valid target to (i.e. a Gaia source ID, a :class:`~ATK.Models.Target`, or an astropy :class:`~astropy.coordinates.SkyCoord`) to :meth:`~ATK.Models.DataPages.open`. For example, to open the above **datapage**:

# sphinx_gallery_start_ignore
doc = Document()
doc.add_root(datapages.figures[0])
# sphinx_gallery_end_ignore
datapages.open(target)
# sphinx_gallery_start_ignore
datapages.figures[0]
# sphinx_gallery_end_ignore

# %%
#
# .. note::
#
#    As with :meth:`~ATK.Models.DataSet.split`, matching a **datapage** to a :class:`~astropy.coordinates.SkyCoord` also requires a ``radius``.
#
# |
# |
#
# Saving Datapages To Local Files
# ===============================
# Saving and Opening
# ------------------
# To open a **datapage** and save it to local files, :meth:`~ATK.Models.DataPages.open` can be provided with a ``path``:
#
# .. code-block:: python
#
#    target = 587316166180416640
#
#    datapages = ...
#
#    datapages.open(target, path="example_datapage.html")
#
# |
# |
#
# Saving Without Opening
# ----------------------
# **Datapages** can also be saved without opening them via :meth:`~ATK.Models.DataPages.save`:
#
# .. code-block:: python
#
#    target = 587316166180416640
#
#    datapages.save(target, path="example_datapage.html")

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
