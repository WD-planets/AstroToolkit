from AstroToolkit.Tools import imagequery,plotimage,plothrd,sedquery,plotsed,spectrumquery,plotspectrum,lightcurvequery,plotlightcurve,plotpowspec,getbuttons,gridsetup,getmdtable,showplot,export
from bokeh.layouts import column,row,layout
from bokeh.plotting import output_file

# source = Hu Leo
#source=2552928187080872832
source=587316166180416640

# set grid size (scales overally size of datapage)
grid_size=275

# get image data and plot it
image_data=imagequery(survey='any',source=source,overlays='gaia,galex_nuv,galex_fuv')
image=plotimage(image_data)

# get hrd
hrd=plothrd(source=source)

# get sed data and plot it
sed_data=sedquery(source=source)
sed=plotsed(sed_data)

# get spectrum data and plot it
spectrum_data=spectrumquery(survey='sdss',source=source)
spectrum=plotspectrum(spectrum_data)

# get lightcurve data [g,r,i]
lightcurve_data=lightcurvequery(survey='ztf',source=source)

# plot each lightcurve in lightcurve_data and set colour of each plot
lightcurves=plotlightcurve(lightcurve_data,colours=['green','red','blue'])

# plot power spectrum data (lightcurve_data)
power_spectrum=plotpowspec(lightcurve_data)

# get SIMBAD and Vizier buttons
buttons=getbuttons(grid_size=grid_size,source=source)

# make a custom metadatatable entry
custom_entry={
    'parameters':['one'],
    'values':['two'],
    'errors':['three'],
    'notes':['four']
    }

# get metadata table
metadata=getmdtable(source=source,metadata={'gaia':'default','galex':'default','panstarrs':'default','skymapper':'default','sdss':'default','wise':'default','twomass':'default','custom':custom_entry})

# get plots, formatted into grid dimensions
grid_plots=gridsetup(dimensions={'width':6,'height':5},plots=[{'name':'image','figure':image,'width':2,'height':2},
                                                              {'name':'hrd','figure':hrd,'width':2,'height':2},
                                                              {'name':'sed','figure':sed,'width':2,'height':1},
                                                              {'name':'lightcurves','figure':lightcurves,'width':2,'height':1},
                                                              {'name':'buttons','figure':buttons,'width':1,'height':1},
                                                              {'name':'spectrum','figure':spectrum,'width':3,'height':1},
                                                              {'name':'powspec','figure':power_spectrum,'width':2,'height':1},
                                                              {'name':'metadata_table','figure':metadata,'width':6,'height':2}],
                                                              grid_size=grid_size)

# set up the final grid
datapage=layout(column(
                row(grid_plots['image'],grid_plots['hrd'],column(grid_plots['sed'],grid_plots['lightcurves'])),
                row(grid_plots['buttons'],grid_plots['spectrum'],grid_plots['powspec']),
                row(grid_plots['metadata_table'])
	            ))

output_file(f'{source}_datapage.html')

showplot(datapage)