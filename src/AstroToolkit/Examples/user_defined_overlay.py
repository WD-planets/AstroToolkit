import os
from pathlib import Path

from AstroToolkit.Tools import query

from .examples_utilities import format, go_to_static

os.chdir(Path(__file__).parent.absolute())

image_data = query(kind="image", source=2552928187080872832, survey="panstarrs", overlays=["allwise"], size=120)

image_data.plot(simbad_search_radius=5)
image_data.plotname = "user_defined_overlay.html"
image_data.figure = format(image_data.figure, 0.5, 0.5)
go_to_static()
image_data.showplot()
