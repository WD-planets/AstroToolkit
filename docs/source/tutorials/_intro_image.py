from pathlib import Path

from _utilities import format_plot
from bokeh.document import Document

from ATK.Tools import query

image_path = Path(__file__).parent / "intro_image.fits.gz"

image = query(kind="image", targets=2552928187080872832, survey="dss1", band="blue", size=120, overlays=["gaia", "galex"], path=image_path)
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
