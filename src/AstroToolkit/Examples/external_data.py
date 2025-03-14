import os
from pathlib import Path

import pandas as pd

from AstroToolkit.Models import CustomLightcurveStruct

data = pd.read_csv(
    os.path.join(Path(__file__).parent.absolute(), "AR_Sco_TNT.txt"), delimiter="\s+"
)

lightcurve = CustomLightcurveStruct(source=6050296829033196032).showdata()

lightcurve.survey = "TNT"
lightcurve.data = [
    {
        "band": "g",
        "hjd": data["mjd"].tolist(),
        "mag": data["flux"].tolist(),
        "mag_err": data["error"].tolist(),
    }
]

lightcurve.plot(colours=["green"]).showplot()
lightcurve.plot(kind="powspec", start_freq=650, stop_freq=800).showplot()
lightcurve.plot(
    kind="phasefold", bins=300, foverlay=False, freq=6.74157303371
).showplot()
