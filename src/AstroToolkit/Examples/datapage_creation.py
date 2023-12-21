from AstroToolkit.Tools import getimage, gethrd,getsed,getsdssspectrum,getztflc,getps,getinfobuttons,getgrid
from bokeh.layouts import column,row,layout
from bokeh.plotting import show,save,output_file

# An example Gaia source_id (*HU Leo)
source=587316166180416640

# Set a grid size to scale all plots to.
grid_size=275

# Grab any plots to include in the data page

# Imaging (uses exhaustitive imaging, i.e. PanSTARRS > SkyMapper > DSS)
image=getimage(source=source)
# HRD
hrd=gethrd(source=source)
# Spectral energy distribution
sed=getsed(source=source)
# SDSS Spectra
spectrum=getsdssspectrum(source=source)
# Lightcurves (return_raw returns a list of lightcurves rather than combining them into one plot)
g_lightcurve,r_lightcurve,i_lightcurve=getztflc(source=source,return_raw=True)
# Power spectrum
power_spectrum=getps(source=source)
# SIMBAD and Vizier buttons
buttons=getinfobuttons(grid_size=grid_size,source=source)

# Defines a grid of size 6 units (of grid_size) long and 3 units (of grid_size) tall, and passes plots into getgrid with their desired unit dimensions.
# The total area covered must equal the area covered by the dimensions (i.e. here we have a total area of 18 units in both the dimensions parameter and the summed plot areas)
plots=getgrid(dimensions=[6,3],plots=[[image,2,2],[buttons,1,1],[g_lightcurve,1,1],[hrd,2,2],[r_lightcurve,1,1],[i_lightcurve,1,1],[sed,2,1],[spectrum,2,1],[power_spectrum,2,1]],grid_size=grid_size)

# getgrid checks for valid input, and then returns these plots stripped of their dimensions. These can then be passed into your Bokeh layout.

# Here, we have a row of 3 columns. The width of each column is set by its widest element.
datapage=layout(
		row(
			column([plots[0],row(plots[1],plots[2])]),
			column([plots[3],row(plots[4],plots[5])]),
			column([plots[6],plots[7],plots[8]])
	))

# Above, the first column consists of a 2x2 panel with a row of two 1x1 panels underneath it. The second column consists of the same, and the third column consists of three 2x1 panels.
# The total size is therefore 6x3 as expected.

# You can then set the name of the output file and show (or save) it.
output_file(f'{source}_datapage.html')
show(datapage)