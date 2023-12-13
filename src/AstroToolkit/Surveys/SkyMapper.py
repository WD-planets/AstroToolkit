import pandas as pd
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import PercentileInterval, AsinhStretch
from io import BytesIO
import PIL
from PIL import Image
import requests
import numpy as np
from astropy.wcs import WCS
import warnings
from astropy.utils.exceptions import AstropyWarning
import math
import astropy
import re
from bokeh.plotting import figure
from bokeh.models import Range1d

#Ignore warning about deprecated header
warnings.simplefilter('ignore', category=AstropyWarning)

def get_image(ra,dec,size=30,band='g'):
	service='https://api.skymapper.nci.org.au/public/siap/dr2/query'
	url=f'{service}?POS={ra},{dec}&SIZE={size},{size}&BAND={band}&FORMAT=image/fits&VERB=3&INTERSECT=covers&RESPONSEFORMAT=CSV'

	try:
		table=pd.read_csv(url)
		table=Table.from_pandas(table)
	except:
		print('[SkyMapper: get_image_data] Error: experiencing issues with SkyMapper.')
		return None, None, None
	
	if len(table)>0:
		urlMain=table['get_image'][0]
	else:
		print('[SkyMapper: get_image_data] Error: no SkyMapper images found using given fields.')
		return None, None, None
	
	# Need to add filter check functionality
	'''
	shortnames=[]
	for i in range(0,len(table)):
		shortnames.append(table['shortname'][i])
	
	missing_filters=[]
	if len(table)!=len(band):
		for i in range(0,len(band)):
			filterCheck='.'+band[i]+'.'
			if any(filterCheck in shortname for shortname in shortnames):
				pass
			else:
				print('[PanSTARRS: get_image_data] Note:  no PanSTARRS data in '+band[i]+' band')
				missing_filters.append(band[i])
				
	for i in range(0,len(missing_filters)):
		band=band.replace(missing_filters[i],'')
	'''
		
	try:
		r=requests.get(urlMain,timeout=15)
		if r.status_code!=200:
			print('[SkyMapper: get_image_data] Error: experiencing issues with SkyMapper.')
			return None, None, None
	except:
		print('[SkyMapper: get_image_data] Error: experiencing issues with SkyMapper.')
		return None, None, None
	
	try:
		fh=fits.open(BytesIO(r.content))
		fitsImage=fh[0].data
		image_header=fh[0].header
		fitsImage[np.isnan(fitsImage)]=0.0
		transform=AsinhStretch()+PercentileInterval(95)
		image=transform(fitsImage)
			
		if not isinstance(image,np.ndarray):
			print('[SkyMapper: get_image_data] Error: no SkyMapper images found using given fields.')
			return None, None, None
			
		return image, image_header, band
	
	except:
		print('[SkyMapper: get_image_data] Error: no SkyMapper images found using given fields.')
		return None, None, None

def get_info(ra,dec,size=30,band='g'):
	ra,dec,size,sizeAS,band=input_management(ra=ra,dec=dec,size=size,band=band)

	image_data=get_image(ra=ra,dec=dec,size=size,band=band)[1]
	
	if not isinstance(image_data,astropy.io.fits.header.Header):
		return None, None
	
	wcs=WCS(image_data)
	
	mjd=image_data['MJD-OBS']
	
	return wcs, mjd

def input_management(ra,dec,size=30,band='g'):
	if len(band)!=len(set(band)):
		raise Exception('duplicate detected in filter selection')

	if not isinstance(ra,float):
		try:
			ra=float(ra)
		except:
			raise Exception('object RA must be a float or an integer')
	if not isinstance(dec,float):
		try:
			dec=float(dec)
		except:
			raise Exception('object DEC must be a float or an integer')
	
	band=band.lower()
	
	if isinstance(size,int):
		pass
	else:
		raise Exception('size must take an integer value')

	if size>600:
		raise Exception('SkyMapper image takes a maximum size of 10 arcmin')
		
	sizeAS=size
	size=size/3600
	
	if re.match('^[grizuv]+$', band):
			pass
	else:
		raise Exception('bands must be in list: [g,r,i,z,u,v]')
	
	band=list(band)
	temp_string=''
	for i in range(0,len(band)):
		temp_string+=(band[i]+',')
	band=temp_string[:-1]

	return ra,dec,size,sizeAS,band

def get_plot(ra,dec,size=30,band='g'):
	ra,dec,size,sizeAS,band=input_management(ra,dec,size=size,band=band)

	# Create axes and plot image
	plot=figure(width=400,height=400,title=f'SkyMapper Image ({sizeAS}")')
	# Fixes scaling issue
	plot.min_border=75

	plot.xgrid.grid_line_color=None
	plot.ygrid.grid_line_color=None	

	image,image_data,_=get_image(ra=ra,dec=dec,size=size,band=band)
	
	if not isinstance(image,np.ndarray) or not isinstance(image_data,astropy.io.fits.header.Header):
		return None

	# Get x and y limit (in pixels) from image header
	xlim=image_data['NAXIS1']
	ylim=image_data['NAXIS2']
	
	# Get points on axes in world coordinates (deg)
	points=[]
	wcs=WCS(image_data)

	# Get array of pixels in x and y
	x_points=np.arange(start=0,stop=xlim+1,step=1)
	y_points=np.arange(start=0,stop=ylim+1,step=1)

	# Apply transformation to coordinates
	coords=wcs.all_pix2world(x_points,y_points,1)
	x_points,y_points=coords[0],coords[1]
	
	# Get image size in degrees
	x_range=max(x_points)-min(x_points)
	y_range=max(y_points)-min(y_points)

	plot.x_range=Range1d(max(x_points),min(x_points))
	plot.y_range=Range1d(min(y_points),max(y_points))

	plot.image(image=[image],x=x_points[0],y=y_points[0],dw=x_range,dh=y_range,palette='Greys256',level="image",origin='bottom_right',anchor='bottom_right')

	return plot

def SkyMapperQueryCoords(ObjRa,ObjDec,radius=3):
	radius=radius/3600

	#Had to use JSON format to get this to work, the CSV format used is not easily formattable by pandas
	url=f'https://skymapper.anu.edu.au/sm-cone/public/query?RA={ObjRa}&DEC={ObjDec}&SR={radius}&RESPONSEFORMAT=CSV'
	try:
		r=requests.get(url)
	except:
		print('[SkyMapper: SkyMapperQueryCoords] Error: experiencing issues with SkyMapper.')
		return None
	
	if r.status_code!=200:
		print('[SkyMapper: SkyMapperQueryCoords] Error: experiencing issues with SkyMapper.')
		return None

	data=BytesIO(r.content)
	try:
		data=pd.read_csv(data)
	except:
		print('[SkyMapper: SkyMapperQueryCoords] Note: No data found using given coordinates')
	
	return data
	
def SkyMapperGetPhotometryCoords(ObjRa,ObjDec,radius=3):
	data=SkyMapperQueryCoords(ObjRa,ObjDec,radius=radius)
	
	if not isinstance(data,pd.DataFrame) or data.empty:
		print('[SkyMapper: SkyMapperGetPhotometryCoords] Note: No SkyMapper data at given coordinates.')
		return None
	
	photometry=data[['raj2000','dej2000','object_id','g_psf','e_g_psf','r_psf','e_r_psf','i_psf','e_i_psf','z_psf','e_z_psf','u_psf','e_u_psf','v_psf','e_v_psf']].copy()
	return photometry