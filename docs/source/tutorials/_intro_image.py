from _utilities import format_plot
from bokeh.document import Document

from ATK.Tools import query

image = query(
    kind="image",
    targets=2552928187080872832,
    survey="dss1",
    band="blue",
    size=120,
    overlays=["gaia", "galex"],
    path="./source/tutorials/intro_image.fits",
)
# sphinx_gallery_start_ignore
image.plot()
figure = format_plot(image.figure, 2, 2)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
