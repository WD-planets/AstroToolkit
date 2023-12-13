import requests
import pandas as pd
from io import BytesIO
import re
import numpy as np
from requests.adapters import HTTPAdapter, Retry
import matplotlib
from bokeh.plotting import figure, show
from bokeh.models import Whisker, ColumnDataSource
from bokeh.transform import linear_cmap
from bokeh.layouts import gridplot

def getData(ra,dec,radius=3):
	#Filters inputs
	if not isinstance(ra,float):
		try:
			ra=float(ra)
		except:
			raise Exception('object RA must be a float or an integer')
	if not isinstance(radius,float):
		try:
			radius=float(radius)
		except:
			raise Exception('search radius must be a float or an integer')
	if not isinstance(dec,float):
		try:
			dec=float(dec)
		except:
			raise Exception('object DEC must be a float or an integer') 
	
	#Convert radius into degrees
	radius=radius/3600
	
	#Perform query
	service='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'
	url=f'{service}?POS=CIRCLE {ra} {dec} {radius}&BANDNAME=g,r,i&FORMAT=CSV'
	
	s=requests.Session()
	retries=Retry(total=5,backoff_factor=0.1,status_forcelist=[500,502,503,504])
	s.mount('http://',HTTPAdapter(max_retries=retries))
	try:
		r=s.get(url,timeout=15)
	except:
		print('[ZTF: getData] Error: experiencing issues with ZTF.')
		return [None, None, None]
	
	if r.status_code!=200:
		print('[ZTF: getData] Error: experiencing issues with ZTF.')
		return [None, None, None]

	#Convert to pandas array (containing all bands)
	data=pd.read_csv(BytesIO(r.content))
	if len(data)==0:
		print('[ZTF: getData] Error:  no ZTF data found using given fields')
		return [None, None, None]
	
	#Split into separate bands
	gData=data.loc[data['filtercode']=='zg']
	rData=data.loc[data['filtercode']=='zr']
	iData=data.loc[data['filtercode']=='zi']
	
	gData,rData,iData=gData.reset_index(drop=True),rData.reset_index(drop=True),iData.reset_index(drop=True)
	
	if gData.empty:
		gData=None
	if rData.empty:
		rData=None
	if iData.empty:
		iData=None

	Data=[gData,rData,iData]
	
	return Data
	
def getLightCurve(ra,dec,radius=3,return_raw=False):
	Data=getData(ra=ra,dec=dec,radius=radius)
	
	try:
		if Data==None:
			if return_raw==False:
				return None
			else:
				return None, None, None
	except:
		pass
	
	#Check if there is any returned data
	if not isinstance(Data[0],pd.DataFrame) and not isinstance(Data[1],pd.DataFrame) and not isinstance(Data[2],pd.DataFrame):
		if return_raw==False:
			return None
		else:
			return None, None, None

	gBand,rBand,iBand=False,False,False
	
	g_colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['greenyellow','forestgreen','greenyellow'])
	g_palette=[matplotlib.colors.rgb2hex(c) for c in g_colourmap(np.linspace(0,1,255))]	

	r_colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['yellow','red','yellow'])
	r_palette=[matplotlib.colors.rgb2hex(c) for c in r_colourmap(np.linspace(0,1,255))]	

	i_colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['aqua','royalblue','aqua'])
	i_palette=[matplotlib.colors.rgb2hex(c) for c in i_colourmap(np.linspace(0,1,255))]	

	#Separate data into different bands, find number of subplots needed and which bands have data
	if isinstance(Data[0],pd.DataFrame):
		gData=Data[0]
		g_hjd=gData.loc[:,'hjd'].tolist()
		g_hjd_min=min(g_hjd)
		g_hjd=[x-g_hjd_min for x in g_hjd]
		g_mag=gData.loc[:,'mag'].tolist()
		g_mag=np.array(g_mag)
		g_mag_err=gData.loc[:,'magerr'].tolist()
		gBand=True
	if isinstance(Data[1],pd.DataFrame):
		rData=Data[1]
		r_hjd=rData.loc[:,'hjd'].tolist()
		r_hjd_min=min(r_hjd)
		r_hjd=[x-r_hjd_min for x in r_hjd]
		r_mag=rData.loc[:,'mag'].tolist()
		r_mag=np.array(r_mag)
		r_mag_err=rData.loc[:,'magerr'].tolist()
		rBand=True
	if isinstance(Data[2],pd.DataFrame):
		iData=Data[2]
		i_hjd=iData.loc[:,'hjd'].tolist()
		i_hjd_min=min(i_hjd)
		i_hjd=[x-i_hjd_min for x in i_hjd]
		i_mag=iData.loc[:,'mag'].tolist()
		i_mag=np.array(i_mag)
		i_mag_err=iData.loc[:,'magerr'].tolist()
		iBand=True
	
	#Plot data
	if gBand:
		gPlot=figure(width=400,height=400,title='ZTF g Lightcurve',x_axis_label=r'\[\text{Observation Date [days]}\]',y_axis_label=r'\[\text{g}\]')
		
		gPlot.y_range.flipped=True

		upper=[x+e for x,e in zip(g_mag,g_mag_err)]
		lower=[x-e for x,e in zip(g_mag,g_mag_err)]
		source=ColumnDataSource(data=dict(hjd=g_hjd,mag=g_mag,upper=upper,lower=lower))

		g_mapper=linear_cmap(field_name='mag',palette=g_palette,low=min(g_mag),high=max(g_mag))

		gPlot.circle(x='hjd',y='mag',source=source,color=g_mapper)
		errors=Whisker(source=source,base='hjd',upper='upper',lower='lower',line_width=0.5,line_color='forestgreen')
		errors.upper_head.line_width,errors.lower_head.line_width=0.5,0.5
		errors.upper_head.size,errors.lower_head.size=3,3
		errors.upper_head.line_color,errors.lower_head.line_color='forestgreen','forestgreen'		

		gPlot.add_layout(errors)
	else:
		gPlot=None

	if rBand:
		rPlot=figure(width=400,height=400,title='ZTF r Lightcurve',x_axis_label=r'\[\text{Observation Date [days]}\]',y_axis_label=r'\[\text{r}\]')
		
		rPlot.y_range.flipped=True

		upper=[x+e for x,e in zip(r_mag,r_mag_err)]
		lower=[x-e for x,e in zip(r_mag,r_mag_err)]
		source=ColumnDataSource(data=dict(hjd=r_hjd,mag=r_mag,upper=upper,lower=lower))

		r_mapper=linear_cmap(field_name='mag',palette=r_palette,low=min(r_mag),high=max(r_mag))

		rPlot.circle(x='hjd',y='mag',source=source,color=r_mapper)
		errors=Whisker(source=source,base='hjd',upper='upper',lower='lower',line_width=0.5,line_color='red')
		errors.upper_head.line_width,errors.lower_head.line_width=0.5,0.5
		errors.upper_head.size,errors.lower_head.size=3,3
		errors.upper_head.line_color,errors.lower_head.line_color='red','red'		

		rPlot.add_layout(errors)
	else:
		rPlot=None
		
	if iBand:
		iPlot=figure(width=400,height=400,title='ZTF i Lightcurve',x_axis_label=r'\[\text{Observation Date [days]}\]',y_axis_label=r'\[\text{i}\]')
		
		iPlot.y_range.flipped=True

		upper=[x+e for x,e in zip(i_mag,i_mag_err)]
		lower=[x-e for x,e in zip(i_mag,i_mag_err)]
		source=ColumnDataSource(data=dict(hjd=i_hjd,mag=i_mag,upper=upper,lower=lower))

		i_mapper=linear_cmap(field_name='mag',palette=i_palette,low=min(i_mag),high=max(i_mag))

		iPlot.circle(x='hjd',y='mag',source=source,color=i_mapper)
		errors=Whisker(source=source,base='hjd',upper='upper',lower='lower',line_width=0.5,line_color='royalblue')
		errors.upper_head.line_width,errors.lower_head.line_width=0.5,0.5
		errors.upper_head.size,errors.lower_head.size=3,3
		errors.upper_head.line_color,errors.lower_head.line_color='royalblue','royalblue'
		
		iPlot.add_layout(errors)
	else:
		iPlot=None

	if return_raw==True:
		return gPlot,rPlot,iPlot
	else:
		grid=gridplot([[gPlot,rPlot,iPlot]],width=400,height=400)
		return grid

if __name__=='__main__':
	plot=getLightCurve(245.4470133276900,-22.8862182469700)
	show(plot)