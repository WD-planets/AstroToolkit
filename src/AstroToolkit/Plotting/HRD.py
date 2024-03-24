from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
import cmasher as cmr
from bokeh.plotting import figure
import os
from importlib_resources import files
from bokeh.models import CustomJS
from bokeh import events

from ..Tools import dataquery

fits_file = files('AstroToolkit.Plotting').joinpath('backdrop_hrd_allmags.fits')

'''
Creates the HRD plot with a source/sources as an overlay
'''
def get_plot(source=None,sources=None):
	# reads the background file ('backdrop_hrd_allmags.fits') into separate datasets based on their object classifications.
	background=Table.read(fits_file).to_pandas()
	background_100pc=background.loc[background['type']==b'100pc']
	background_wd=background.loc[background['type']==b'WD']
	background_cv=background.loc[background['type']==b'CV']
	background_wd_dm=background.loc[background['type']==b'WD+dM']
	background_sd=background.loc[background['type']==b'SD']
	
	# combines these into a split dataset
	background_split=[background_100pc,background_wd,background_wd_dm,background_cv,background_sd]
	
	x_arr=[]
	y_arr=[]
	colour_arr=[]
	
	# calculates x (bp-rp) and y (abs g) positions of all objects in the background sample
	for sample in background_split:
		x=sample['phot_bp_mean_mag'].values-sample['phot_rp_mean_mag'].values
		y=sample['phot_g_mean_mag'].values+5*np.log10(sample['parallax'].values/1000)+5
		
		x_arr.append(x)
		y_arr.append(y)
	
	# Sets the colour maps used by each category of object
	colour_maps=['Greens','Purples','Reds','Blues','YlOrBr']
	
	# Removes extreme ends of colour maps for visibility (i.e. strips colours that are too light)
	for i in range(0,len(colour_maps)):
		colour_maps[i]=cmr.get_sub_cmap(colour_maps[i],0.4,1.0)

	# Creates the effect of darker colours in higher density regions of each category
	for i in range(0,len(x_arr)):
		#stacks x and y values to calculate a gaussian density and then uses this to normalise across the colour map based on density of nearby data points
		values=np.vstack([x_arr[i],y_arr[i]])
		kernel=stats.gaussian_kde(values)
		weights=kernel(values)
		weights=weights/weights.max()
		
		# uses a hex representation of RGB colour values to create the final colour maps weighted for each object based on nearby point density
		colour=["#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b, _ in 255*colour_maps[i](weights)]
		colour_arr.append(colour)

	# Create plot
	plot=figure(width=400,height=400,title='Gaia HRD',x_axis_label=r'\[\text{bp-rp}\]',y_axis_label=r'\[\text{abs g}\]')
	plot.y_range.flipped = True

	# sets labels for legend
	labels=['100pc','WDs','CVs','WD+dM','SDs']
	# places background points with their colours set, and addstheir respective category to the legend
	for i in range(0,len(x_arr)):
		plot.scatter(x_arr[i],y_arr[i],size=3,color=colour_arr[i],legend_label=labels[i])

	# overlays source on top of background sample as long as all the necessary data is available
	if source!=None:
		data=dataquery(survey='gaia',source=source)['data']
		if data==None:
			return None
		gmag,bpmag,rpmag,parallax=data['phot_g_mean_mag'][0],data['phot_bp_mean_mag'][0],data['phot_rp_mean_mag'][0],data['parallax'][0]
		x=bpmag-rpmag
		y=gmag+5*np.log10(parallax/1000)+5
		
		if not math.isnan(gmag) and not math.isnan(bpmag) and not math.isnan(rpmag) and not math.isnan(parallax):
			plot.scatter(x=x,y=y,marker='square_dot',line_color='black',fill_color=None,size=20,line_width=2,legend_label='Source')
		else:
			return None
	
	# same as above but handles multiple sources
	elif sources!=None:
		for i in range(0,len(sources)):
			data=dataquery(survey='gaia',source=sources[i])['data']
			gmag,bpmag,rpmag,parallax=data['phot_g_mean_mag'][0],data['phot_bp_mean_mag'][0],data['phot_rp_mean_mag'][0],data['parallax'][0]
			x=bpmag-rpmag
			y=gmag+5*np.log10(parallax/1000)+5
			
			if not math.isnan(gmag) and not math.isnan(bpmag) and not math.isnan(rpmag) and not math.isnan(parallax):
				plot.scatter(x=x,y=y,marker='square_dot',line_color='black',fill_color=None,size=20,line_width=2,legend_label='Source')
			else:
				pass

	if plot!=None:
		plot.legend.click_policy="hide"		

		# Double click to hide legend
		toggle_legend_js = CustomJS(args=dict(leg=plot.legend[0]), code='''
			 if (leg.visible) {
				 leg.visible = false
				 }
			 else {
				 leg.visible = true
			 }
		''')
	
		plot.js_on_event(events.DoubleTap, toggle_legend_js) 

	return plot