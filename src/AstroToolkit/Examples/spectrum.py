from AstroToolkit.Tools import query

# specify a Gaia source
source = 587316166180416640

# retrieve an SDSS spectrum for the gaia source, plot it and then show it.
# As no radius is given, this will be taken from the config.
spectrum = query(kind="spectrum", source=source, survey="sdss").plot().showplot()
