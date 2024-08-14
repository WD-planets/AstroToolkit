from AstroToolkit.Tools import query

# specify a Gaia source
source = 587316166180416640

# retrieve ZTF light curve data and plot it, specifying the colours for each band.
# No radius  is given, so this will be taken from the config.
figure = query(kind="lightcurve", source=source, survey="ztf").plot(
    colours=["green", "red", "blue"]
)

# show the lightcurves in the browser (and save to a static .html file)
figure.showplot()
