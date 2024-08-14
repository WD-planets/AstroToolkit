from AstroToolkit.Tools import query

# specify a Gaia source
source = 587316166180416640

# retrieve and plot an sed, with an overlayed SDSS spectrum.
sed = query(kind="sed", source=source).plot(spectrum_overlay=True, survey="sdss")

# show the figure in the browser (and save to a static .html file)
sed.showplot()
