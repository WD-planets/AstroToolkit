import os
from pathlib import Path

from AstroToolkit.Tools import query, readdata

os.chdir(Path(__file__).parent.absolute())

gaia_data = query(kind="data", source=587316166180416640, survey="galex")
gaia_data.savedata("test_data.fits")
data = readdata("test_data.fits")
data.showdata()

lightcurve_data = query(
    kind="lightcurve",
    source=6050296829033196032,
    survey="ztf",
    check_exists="test_lightcurve.fits",
)
