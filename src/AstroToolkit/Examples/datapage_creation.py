import os
from pathlib import Path

from AstroToolkit.Datapages import buttons, datapage, datatable
from AstroToolkit.Tools import query

from .examples_utilities import go_to_static

file_dir = Path(__file__).parent.absolute()
os.chdir(Path(__file__).parent.absolute())

source = 587316166180416640

image = query(
    kind="image", survey="panstarrs", source=source, overlays=["gaia", "galex"], check_exists="datapage_image"
).plot()

hrd = query(kind="hrd", sources=source, check_exists=os.path.join(file_dir, "datapage_hrd")).plot()

spectrum = query(
    kind="spectrum", survey="sdss", source=source, check_exists=os.path.join(file_dir, "datapage_spectrum")
).plot()

sed = query(kind="sed", source=source, check_exists="datapage_sed").plot(spectrum_overlay=spectrum)

lightcurves = query(
    kind="lightcurve", survey="ztf", source=source, check_exists=os.path.join(file_dir, "datapage_lightcurve")
).plot(colours=["green", "red", "blue"])

powspec = query(
    kind="lightcurve", survey="ztf", source=source, check_exists=os.path.join(file_dir, "datapage_lightcurve")
).plot(kind="powspec")

buttons = buttons(source=source)

metadata = datatable(
    source=source,
    entries=[
        {"kind": "atk_defaults", "surveys": ["gaia", "galex", "panstarrs", "skymapper", "sdss", "wise", "twomass"]}
    ],
)

datapage = datapage(
    dimensions={"width": 6, "height": 6},
    panels=[
        {"name": "image", "figure": image, "width": 2, "height": 2},
        {"name": "hrd", "figure": hrd, "width": 2, "height": 2},
        {"name": "sed", "figure": sed, "width": 3, "height": 2},
        {"name": "buttons", "figure": buttons, "width": 2, "height": 2},
        {"name": "lightcurves", "figure": lightcurves, "width": 4, "height": 2},
        {"name": "powspec", "figure": powspec, "width": 3, "height": 2},
        {"name": "spectrum", "figure": spectrum, "width": 5, "height": 2},
        {"name": "metadata_table", "figure": metadata, "width": 7, "height": 2},
    ],
    layout=[["image", "sed", "buttons"], ["hrd", "spectrum"], ["lightcurves", "powspec"], ["metadata_table"]],
)

go_to_static()
datapage.showplot(f"{source}_datapage")
