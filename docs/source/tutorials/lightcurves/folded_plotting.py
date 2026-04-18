"""
################################
Folded Light Curve Customisation
################################
When plotting *folded* light curves, some additional arguments are available.

|

"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query

asassn_query = query("lightcurve", targets=5346631922949364864, survey="asassn", path="example_lightcurve.fits.gz")
folded_data = asassn_query.apply("fold", freq=2.773227732277323, inplace=False)

# %%
# Relative Magnitudes
# ===================
# Bands can be magnitude-subtracted (and therefore aligned) by passing ``subtract = 'median'`` (default) or ``subtract = 'mean'``.

# sphinx_gallery_start_ignore
folded_data.plot(subtract="mean")
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open(subtract="mean")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# To disable this behavour, pass ``subtract = None``.

# sphinx_gallery_start_ignore
folded_data.plot(subtract=None)
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open(subtract=None)
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
#
# Aligning in Phase
# =================
# Folded light curves can be shifted in phase to align with a chosen feature. By default, light curves are assinged to a local maximum (i.e. ``align = 'max'``). To instead align with a local minimum, pass ``align = 'min'``:

# sphinx_gallery_start_ignore
folded_data.plot(align="min")
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open(align="min")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# Alternatively, pass ``align = 'median'`` or ``align = 'mean'`` to align to the median or mean magnitude, respectively:

# sphinx_gallery_start_ignore
folded_data.plot(align="median", subtract=None)
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open(align="median")
# sphinx_gallery_start_ignore
figure
# sphinx_gallery_end_ignore

# %%
# |
# 
# Setting Repetitions
# ===================
# By default, folded light curves are repeated to show two complete modulations in brightness. This can be changed by passing ``repeat`` with a number of repetitions:

# sphinx_gallery_start_ignore
folded_data.plot(repeat=1)
figure = format_plot(folded_data.figure, 3, 1.5)
doc = Document()
doc.add_root(figure)
# sphinx_gallery_end_ignore
folded_data.open(repeat=1)
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
