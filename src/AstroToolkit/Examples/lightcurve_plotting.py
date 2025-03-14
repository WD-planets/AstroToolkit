import os
from pathlib import Path

from AstroToolkit.Tools import query

from .examples_utilities import format, go_to_static

os.chdir(Path(__file__).parent.absolute())

lightcurve_data = query(
    kind="lightcurve",
    source=6050296829033196032,
    survey="ztf",
    check_exists="lightcurve.fits",
)

go_to_static()

lightcurve_data.showdata()
lightcurve_data.plot()
lightcurve_data.figure = format(lightcurve_data.figure, change_size=False)
lightcurve_data.plotname = "lightcurve1.html"
lightcurve_data.showplot()

lightcurve_data.plot(bands=["g", "r", "i"], colours=["green", "red", "blue"])
lightcurve_data.figure = format(lightcurve_data.figure, change_size=False)
lightcurve_data.plotname = "lightcurve2.html"
lightcurve_data.showplot()

lightcurve_data.plot(kind="powspec")
lightcurve_data.figure = format(lightcurve_data.figure, 0.5, 0.5)
lightcurve_data.plotname = "powspec.html"
lightcurve_data.showplot()

lightcurve_data.plot(kind="phasefold")
lightcurve_data.figure = format(lightcurve_data.figure, 0.5, 0.5)
lightcurve_data.plotname = "phasefold1.html"
lightcurve_data.showplot()

lightcurve_data.plot(kind="phasefold", bins=100, shift=0.115)
lightcurve_data.figure = format(lightcurve_data.figure, 0.5, 0.5)
lightcurve_data.plotname = "phasefold2.html"
lightcurve_data.showplot()
