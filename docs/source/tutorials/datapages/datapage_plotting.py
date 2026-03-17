"""
#########
DataPages
#########
ATK supports the creation of :class:`DataPages <ATK.structures.DataPage.DataPage>` as a means of neatly combining plots into a single grid. To demonstrate this, some data must first be acquired:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query


target = 587316166180416640

lc = query("lightcurve", survey="asassn", targets=target, path="datapage_lc.fits")
pspec = lc.apply("pspec", min=0, max=60, samples=100000, inplace=False)
fold = lc.apply("fold", min=0, max=60, samples=100000, inplace=False)
image = query("image", survey="panstarrs", band="g", targets=target, path="datapage_image.fits")
spec = query("spectrum", survey="sdss", targets=target, path="datapage_spec.fits")
sed = query("sed", targets=target, path="datapage_sed.fits")
hrd = query("hrd", targets=target, path="datapage_hrd.fits")
table = query("datatable", rows={"gaia": ["Gmag", "BPmag", "RPmag"]}, targets=target, path="datapage_table.fits")

# %%
# The only other requirement is to define the :class:`~ATK.structures.DataPage.DataPage` layout. The easiest way to do this is by passing a list of rows of :class:`DataSets <ATK.Models.DataSet>`:

layout = [[image, image, hrd,   hrd,   lc,    lc,    lc,    lc   ],
          [image, image, hrd,   hrd,   lc,    lc,    lc,    lc   ],
          [pspec, pspec, pspec, pspec, fold,  fold,  fold,  fold ],
          [pspec, pspec, pspec, pspec, fold,  fold,  fold,  fold ],
          [spec,  spec,  spec,  spec,  sed,   sed,   sed,   sed  ],
          [spec,  spec,  spec,  spec,  sed,   sed,   sed,   sed  ],
          [table, table, table, table, table, table, table, table],
          [table, table, table, table, table, table, table, table]]

# %%
# This defines a datapage where the first row contains a 2x2 image in the top-left corner, followed by a 2x2 HRD and a 4x2 light curve, etc. Passing this to :func:`~ATK.Visualisation.grid` returns a :class:`~ATK.structures.DataPage.DataPage`:

from ATK import grid

datapage = grid(layout)
datapage.show()

# %%
# The only remaining step is to open or save our :class:`~ATK.structures.DataPage.DataPage`. Since :class:`DataPages <ATK.structures.DataPage.DataPage>` are designed to work with multiple targets, we can open the :class:`~ATK.structures.DataPage.DataPage` for our specific target via one of the :class:`~ATK.structures.DataPage.DataPage`'s open method (these work in much the same way as those of :class:`DataSets <ATK.Models.DataSet>`). Because a Gaia Source ID was used for targeting, we can use the :meth:`~ATK.structures.DataPage.DataPage.open_by_id` method of the returned :class:`~ATK.structures.DataPage.DataPage`:

# sphinx_gallery_start_ignore
doc = Document()
doc.add_root(datapage.figures[0])
# sphinx_gallery_end_ignore
datapage.open_by_id(target)
# sphinx_gallery_start_ignore
datapage.figures[0]
# sphinx_gallery_end_ignore
