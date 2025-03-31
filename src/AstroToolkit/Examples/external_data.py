import os
from pathlib import Path

import pandas as pd

from AstroToolkit.Models import CustomLightcurveStruct

from .examples_utilities import format, go_to_static

os.chdir(Path(__file__).parent.absolute())

data = pd.read_csv(os.path.join(Path(__file__).parent.absolute(), "AR_Sco_TNT.txt"), delimiter="\s+")

lightcurve = CustomLightcurveStruct(source=6050296829033196032).showdata()

lightcurve.survey = "TNT"
lightcurve.data = [
    {"band": "g", "mjd": data["mjd"].tolist(), "flux": data["flux"].tolist(), "flux_err": data["error"].tolist()}
]

go_to_static()

lightcurve.plot(colours=["green"])
lightcurve.plotname = "AR_Sco_TNT_Lightcurve.html"
lightcurve.figure = format(lightcurve.figure, change_size=False)
lightcurve.showplot()

lightcurve.plot(kind="powspec", start_freq=650, stop_freq=800)
lightcurve.plotname = "AR_Sco_TNT_Powspec.html"
lightcurve.figure = format(lightcurve.figure, 0.5, 0.5)
lightcurve.figure.legend.visible = False
lightcurve.figure.xaxis.ticker.desired_num_ticks = 4
lightcurve.showplot()

lightcurve.plot(kind="phasefold", bins=300, foverlay=False, freq=6.74157303371)
lightcurve.plotname = "AR_Sco_TNT_Phasefold.html"
lightcurve.figure = format(lightcurve.figure, change_size=False)
lightcurve.figure.width = 1000
lightcurve.figure.height = 500
lightcurve.showplot()
