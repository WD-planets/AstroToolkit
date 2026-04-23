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
# The only remaining step is to open or save the **datapage**. Like everything in ATK, :class:`DataPages <ATK.Models.DataPages>` are designed to work with multiple targets. A **datapage** for a specific target can therefore be openened via one of the :class:`~ATK.Models.DataPages` object's open methods.
#
# Opening a Datapage by ID 
# ----------------------------
# Since a Gaia Source ID was used for targeting, the :meth:`~ATK.Models.DataPages.open_by_id` method of the returned :class:`~ATK.Models.DataPages` can be used:

# sphinx_gallery_start_ignore
doc = Document()
doc.add_root(datapages.figures[0])
# sphinx_gallery_end_ignore
datapages.open(target)
# sphinx_gallery_start_ignore
datapages.figures[0]
# sphinx_gallery_end_ignore

# %%
# |
#
# Opening a Datapage by Coordinates
# ---------------------------------
# :meth:`~ATK.Models.DataPages.open_by_coord` is an equivalent methods for fetching by coordinates. This opens any **datapages** that target stars within a given ``radius``:
#
# .. code-block:: python
#
#    from astropy.coordinates import SkyCoord
#    import astropy.units as u
#
#    coord = SkyCoord(ra=141.185, dec=8.031, unit="deg", frame="icrs")
#
#    datapages = ...
#
#    datapages.open_by_coord(coord, radius=3*u.arcsec)

# %%
# |
#
# Opening a Datapage by Target
# ----------------------------
# Finally, :meth:`~ATK.Models.DataPages.open_by_target` can be used to open a **datapage** using a :class:`~ATK.Models.Target` directly. This allows for exact matching without the need of a ``radius``:
#
# .. code-block:: python
#
#    from astropy.coordinates import SkyCoord
#
#    from ATK.Models import Target
#
#    coord = SkyCoord(ra=141.185, dec=8.031, unit="deg", frame="icrs")
#    target = Target.from_coord(coord)
#
#    datapages = ...
#
#    datapages.open_by_target(target)
#
# |
# |
#
# Saving Datapages To Local Files
# ===============================
# Saving and Opening
# ------------------
# To open a **datapage** and save it to local files, any of the above methods can be provided with a ``path`` - e.g.
#
# .. code-block:: python
#
#    target = 587316166180416640
#
#    datapages = ...
#
#    datapages.open_by_id(target, path="example_datapage.html")
#
# |
# |
#
# Saving Without Opening
# ----------------------
# **Datapages** can also be saved without opening them via any of the equivalent save methods.
#
# By ID:
#
# .. code-block:: python
#
#    target = 587316166180416640
#
#    datapages.save_by_id(target, path="example_datapage.html")
#
# |
#
# By coordinate:
#
# .. code-block:: python
#
#    from astropy.coordinates import SkyCoord
#
#    coord = SkyCoord(ra=141.185, dec=8.031, unit="deg", frame="icrs")
#
#    datapages = ...
#
#    datapages.save_by_coord(coord, path="example_datapage.html")
#
# |
#
# By :class:`~ATK.Models.Target`:
#
# .. code-block:: python
#
#    from astropy.coordinates import SkyCoord
#
#    from ATK.Models import Target
#
#    coord = SkyCoord(ra=141.185, dec=8.031, unit="deg", frame="icrs")
#    target = Target.from_coord(coord)
#
#    datapages = ...
#
#    datapages.open_by_target(target, path="example_datapage.html")

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
