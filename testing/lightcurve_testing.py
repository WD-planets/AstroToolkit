import ATK.Tools as ATK

SOURCE = 587316166180416640

data = ATK.query("lightcurve", survey="ztf", target=SOURCE)
# data.show()
