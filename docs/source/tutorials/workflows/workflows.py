"""
Putting Everything Together
===========================
Putting :doc:`multi-target operations <../multiple_targets/multi_target_query>`, :meth:`~ATK.Models.DataSet.split`, :meth:`~ATK.Models.DataSet.merge`, and :doc:`datapage generation <../datapages/datapage_plotting>` together allows for sophisticated workflows.

The following example shows a set of 2-target light curve queries to ASAS-SN and TESS:
"""

# sphinx_gallery_start_ignore
# fmt: off
# isort: skip_file
from _utilities import format_plot
from bokeh.document import Document
# sphinx_gallery_end_ignore
from ATK import query, grid

import astropy.units as u

targets = [5346631922949364864, 587316166180416640]

asassn_data = query("lightcurve", targets=targets, survey="asassn", path="assasn_data.fits.gz")
tess_data = query("lightcurve", targets=targets, survey="tess", path="tess_data.fits.gz")

# sphinx_gallery_start_ignore
pass
# sphinx_gallery_end_ignore

# %%
# |
#
# The resulting :class:`DataSets <ATK.Models.DataSet>` can then be combined into a single :class:`~ATK.Models.DataSet`:

all_lcs = asassn_data.merge(tess_data)
all_lcs.show()

# %%
# |
#
# Power spectra can be generated across the merged data set, with all survey/band/target management being performed internally:

all_pspec = all_lcs.apply("pspec", fmin=0, fmax=10, samples=50000, inplace=False)
all_pspec.show()

# %%
# |
#
# The light curves in the merged :class:`~ATK.Models.DataSet` can also be phase-folded and binned:

all_fold = all_lcs.apply("fold", fmin=0, fmax=10, samples=50000, inplace=False)
all_fold.apply("bin", bins=200)
all_fold.show()

# %%
# |
# 
# Finally, the acquired light curve, power spectrum and phase-folded light curve :class:`DataSets <ATK.Models.DataSet>` can be combined into **datapages**:

layout = [[all_lcs,   all_lcs,   all_lcs,  all_lcs,  all_lcs],
          [all_lcs,   all_lcs,   all_lcs,  all_lcs,  all_lcs],
          [all_pspec, all_pspec, all_fold, all_fold, all_fold],
          [all_pspec, all_pspec, all_fold, all_fold, all_fold]]

datapages = grid(layout)

# sphinx_gallery_start_ignore
doc = Document()
plot_map = datapages._plot_map
from ATK.Models import Target
target = Target.from_id(5346631922949364864)
plot_id = plot_map[target._key]
plot = [fig for fig in datapages.figures if fig.id is plot_id][0]
doc.add_root(plot)
# sphinx_gallery_end_ignore
datapages.open_by_id(targets[0])
# sphinx_gallery_start_ignore
plot
# sphinx_gallery_end_ignore

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
