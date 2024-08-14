from AstroToolkit.Tools import query, readdata

# specify a Gaia source
source = 587316166180416640

# retrieve ZTF light curve data for our source and save this structure to a local file
filename = query(kind="lightcurve", source=source, survey="ztf").savedata()

# recreate the original data structure from the local file
recreated_data = readdata(filename)

# plot only the g band of this data in the colour green, and show it.
recreated_data.plot(colours=["green"], bands=["g"]).showplot()
