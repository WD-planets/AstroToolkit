from AstroToolkit.Tools import query

# specify a Gaia source
source1 = 587316166180416640

# Specify a second Gaia source. This is not necessary - any number of sources can be overlayed.
source2 = 6050296829033196032

# retrieve Gaia hrd data for our list of sources, plot it and then show it.
hrd = query(kind="hrd", sources=[source1, source2]).plot().showplot()
