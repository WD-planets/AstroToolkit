from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
import cmasher as cmr
from bokeh.plotting import figure,show
import os
from importlib_resources import files

from ..Surveys.Gaia import GaiaQueryDesignation

fits_file = files('AstroToolkit.Figures').joinpath('backdrop_hrd_allmags.fits')

def get_plot(source=None,sources=None):
	background=Table.read(fits_file).to_pandas()
	background_100pc=background.loc[background['type']==b'100pc']
	background_wd=background.loc[background['type']==b'WD']
	background_cv=background.loc[background['type']==b'CV']
	background_wd_dm=background.loc[background['type']==b'WD+dM']
	background_sd=background.loc[background['type']==b'SD']
	
	background_split=[background_100pc,background_wd,background_wd_dm,background_cv,background_sd]
	
	x_arr=[]
	y_arr=[]
	colour_arr=[]
	
	for sample in background_split:
		x=sample['phot_bp_mean_mag'].values-sample['phot_rp_mean_mag'].values
		y=sample['phot_g_mean_mag'].values+5*np.log10(sample['parallax'].values/1000)+5
		
		x_arr.append(x)
		y_arr.append(y)
	
	colour_maps=['Greens','Purples','Reds','Blues','YlOrBr']
	
	for i in range(0,len(colour_maps)):
		colour_maps[i]=cmr.get_sub_cmap(colour_maps[i],0.4,1.0)

	for i in range(0,len(x_arr)):
		values=np.vstack([x_arr[i],y_arr[i]])
		kernel=stats.gaussian_kde(values)
		weights=kernel(values)
		weights=weights/weights.max()
		
		colour=["#%02x%02x%02x" % (int(r), int(g), int(b)) for r, g, b, _ in 255*colour_maps[i](weights)]
		colour_arr.append(colour)

	plot=figure(width=400,height=400,title='Gaia HRD',x_axis_label=r'\[\text{bp-rp}\]',y_axis_label=r'\[\text{abs g}\]')
	plot.y_range.flipped = True

	labels=['100pc','WDs','CVs','WD+dM','SDs']
	for i in range(0,len(x_arr)):
		plot.scatter(x_arr[i],y_arr[i],size=3,color=colour_arr[i],legend_label=labels[i])

	if source!=None:
		data=GaiaQueryDesignation(Designation=source)
		gmag,bpmag,rpmag,parallax=data['phot_g_mean_mag'].values,data['phot_bp_mean_mag'].values,data['phot_rp_mean_mag'].values,data['parallax'].values
		x=bpmag-rpmag
		y=gmag+5*np.log10(parallax/1000)+5
		
		if not math.isnan(gmag) and not math.isnan(bpmag) and not math.isnan(rpmag) and not math.isnan(parallax):
			plot.scatter(x=x,y=y,marker='square_dot',line_color='black',fill_color=None,size=20,line_width=2,legend_label='Source')
		else:
			return None

	elif sources!=None:
		for i in range(0,len(sources)):
			data=GaiaQueryDesignation(Designation=sources[i])
			gmag,bpmag,rpmag,parallax=data['phot_g_mean_mag'].values,data['phot_bp_mean_mag'].values,data['phot_rp_mean_mag'].values,data['parallax'].values
			x=bpmag-rpmag
			y=gmag+5*np.log10(parallax/1000)+5
			
			if not math.isnan(gmag) and not math.isnan(bpmag) and not math.isnan(rpmag) and not math.isnan(parallax):
				plot.scatter(x=x,y=y,marker='square_dot',line_color='black',fill_color=None,size=20,line_width=2,legend_label='Source')
			else:
				pass

	plot.legend.click_policy="hide"

	return plot