from AstroToolkit.Tools import query

lightcurve_data = query(kind="lightcurve", source=6050296829033196032, survey="ztf")
lightcurve_data.showdata()
lightcurve_data.showplot()

lightcurve_data.plot(bands=["g", "r", "i"], colours=["green", "red", "blue"])
lightcurve_data.showplot()
lightcurve_data.plot(kind="powspec").showplot()
lightcurve_data.plot(kind="phasefold").showplot()
