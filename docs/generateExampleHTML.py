import os
from pathlib import Path

from bokeh.io import output_file, save
from bokeh.models import PanTool

from AstroToolkit.Tools import query

base_width = 700
base_height = 700
base_path = os.path.join(Path(__file__).parent.absolute(), "_static")


def format(plot, width=None, height=None, change_size=True):
    if change_size:
        plot.width = int(width * base_width)
        plot.height = int(height * base_height)
        plot.frame_width = int(width * base_width)
        plot.frame_height = int(height * base_height)
    plot.toolbar_location = None
    plot.grid.grid_line_color = None
    plot.min_border = 0
    plot.title = None
    for tool in plot.select(PanTool):
        plot.remove_tools(tool)
    return plot


data = query(
    kind="image",
    source=2552928187080872832,
    survey="panstarrs",
    overlays=["gaia", "galex"],
    size=120,
    check_exists="image_overlays.fits",
)
data.plot(simbad_search_radius=5)
data.figure = format(data.figure, 0.5, 0.5)
output_file(os.path.join(base_path, "image.html"))
save(data.figure)

data = query(
    kind="lightcurve",
    source=6050296829033196032,
    survey="ztf",
    radius=3,
    check_exists="lightcurve.fits",
)
data.plot()
data.figure = format(data.figure, change_size=False)
output_file(os.path.join(base_path, "lightcurve1.html"))
save(data.figure)

data.plot(bands=["g", "r", "i"], colours=["green", "red", "blue"])
data.figure = format(data.figure, change_size=False)
output_file(os.path.join(base_path, "lightcurve2.html"))
save(data.figure)

data.plot(kind="powspec")
data.figure = format(data.figure, 0.5, 0.5)
output_file(os.path.join(base_path, "powspec.html"))
save(data.figure)

data.plot(kind="phasefold")
data.figure = format(data.figure, 0.5, 0.5)
output_file(os.path.join(base_path, "phasefold1.html"))
save(data.figure)

data.plot(kind="phasefold", bins=100, shift=0.115)
data.figure = format(data.figure, 0.5, 0.5)
output_file(os.path.join(base_path, "phasefold2.html"))
save(data.figure)

data = query(
    kind="lightcurve",
    pos=[141.185, 8.031],
    survey="ztf",
    check_exists="commandline_lightcurve.fits",
).plot()
data.figure = format(data.figure, change_size=False)
output_file(os.path.join(base_path, "commandline_lightcurve.html"))
save(data.figure)

# AR Sco TNT Example

import pandas as pd

from AstroToolkit.Models import CustomLightcurveStruct

data = pd.read_csv(
    "../src/AstroToolkit/Examples/AR_Sco_TNT.txt", delimiter="\s+", dtype=float
)

lightcurve = CustomLightcurveStruct(source=6050296829033196032)

lightcurve.data = [
    {
        "band": "g",
        "hjd": data["mjd"].tolist(),
        "mag": data["flux"].tolist(),
        "mag_err": data["error"].tolist(),
    }
]
lightcurve.survey = "TNT"

lightcurve.plot(colours=["green"])
lightcurve.figure = format(lightcurve.figure, change_size=False)
output_file(os.path.join(base_path, "AR_Sco_TNT_Lightcurve.html"))
save(lightcurve.figure)

lightcurve.plot(kind="powspec", start_freq=650, stop_freq=800)
lightcurve.figure = format(lightcurve.figure, 0.5, 0.5)
lightcurve.figure.xaxis.ticker.desired_num_ticks = 3
output_file(os.path.join(base_path, "AR_Sco_TNT_Powspec.html"))
save(lightcurve.figure)

lightcurve.plot(kind="phasefold", bins=300, foverlay=False, freq=6.74157303371)
lightcurve.figure = format(lightcurve.figure, change_size=False)
figure = lightcurve.figure
figure.width = 1000
figure.height = 500
output_file(os.path.join(base_path, "AR_Sco_TNT_Phasefold.html"))
save(figure)
