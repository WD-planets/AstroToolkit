from AstroToolkit.Tools import query, tsanalysis

# perform PyAOV time series analysis on a light curve data structure
tsanalysis(query(kind="lightcurve", source=6050296829033196032, survey="ztf"))
