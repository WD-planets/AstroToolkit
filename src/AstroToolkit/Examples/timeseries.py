from AstroToolkit.Tools import query

# specify a Gaia source
source = 6050296829033196032

# retrieve ztf data for our Gaia source, plot it as a power spectrum, and then show it.
power_spectrum = (
    query(kind="lightcurve", source=source, survey="ztf")
    .plot(kind="powspec")
    .showplot()
)

# retrieve ztf data for our Gaia source, phase fold it, and then show it.
phase_fold = (
    query(kind="lightcurve", source=source, survey="ztf").plot(kind="phase").showplot()
)
