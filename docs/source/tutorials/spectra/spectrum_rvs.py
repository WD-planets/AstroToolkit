"""
#############################
Detecting Radial Velocities *
#############################
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

sdss_query = query("spectrum", targets=587316166180416640, survey="sdss", path="example_spectrum.fits")
sdss_query.apply("rv_fit", smoothing=5)
sdss_query.show(show_types=True)
# sphinx_gallery_start_ignore
figure = format_plot(sdss_query.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
sdss_query.open()
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
