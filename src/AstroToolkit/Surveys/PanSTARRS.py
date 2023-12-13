from astropy.table import Table
from astropy.io import fits
import astropy
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.wcs import WCS
import PIL
from PIL import Image
from io import BytesIO
import numpy as np
import requests
import re
import warnings
from astropy.utils.exceptions import AstropyWarning
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import Range1d

#Ignore warning about deprecated header
warnings.simplefilter('ignore', category=AstropyWarning)

def get_image(ra,dec,size=30,band='g'):
	# Construct URL for fetching a list of image urls
	service='https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
	url=f'{service}?ra={ra}&dec={dec}&band={band}'
	
	# Try to read the resulting table
	try:
		table=Table.read(url, format='ascii')
	except:
		print('[PanSTARRS: get_image_data] Error: experiencing issues with Pan-STARRS.')
		return None, None, None
	
	if not len(table)>0:
		print('[PanSTARRS: get_image_data] Error: no PanSTARRS images found using given fields.')
		return None, None, None
	
	# Finds any missing filters in returned data and removes them from the 'band' string to signify this
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

	# Constructs the URL used to fetch the images
	url=(f'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={ra}&dec={dec}&size={size}&format=fits')
	
	filterList=['yzirg'.find(x) for x in table['filter']]
	table=table[np.argsort(filterList)]
	
	urlbase=url+'&red='
	url=[]
	for filename in table['filename']:
		url.append(urlbase+filename)

	# Grab images using constructed URL
	try:
		r=requests.get(url[0],timeout=15)
		if r.status_code!=200:
			print('[PanSTARRS: getFitsImage] Error:  experiencing issues with Pan-STARRS.')
			return None, None, None
	except:
		print('[PanSTARRS: get_image_data] Error: experiencing issues with Pan-STARRS.')
		return None, None, None
	
	# Read and format image data
	try:
		fh=fits.open(BytesIO(r.content))
		fitsImage=fh[0].data
		image_header=fh[0].header
		fitsImage[np.isnan(fitsImage)]=0.0
		transform=AsinhStretch()+PercentileInterval(95)
		image=transform(fitsImage)

		return image, image_header, band
	except:
		print('[PanSTARRS: get_image_data] Error: no PanSTARRS images found using given fields. 2')
		return None, None, None

'''
Function for returning the wcs (transformation) and mjd of an image
'''
def get_info(ra,dec,size=30,band='g'):
	# Verify inputs
	ra,dec,size,sizeAS,band=input_management(ra=ra,dec=dec,size=size,band=band)
	# Get image header
	image_data=get_image(ra=ra,dec=dec,size=size,band=band)[1]
	
	if not isinstance(image_data,astropy.io.fits.header.Header):
		return None, None
	
	wcs=WCS(image_data)
	mjd=image_data['MJD-OBS']
	
	return wcs, mjd

'''
Verifies all inputs and transforms them into the units required by queries
'''
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
	
	if not isinstance(size,int):
		try:
			size=int(size)
		except:
			raise Exception('size must take an integer value')
	
	# Adjusts size to unit required by PanSTARRS (0.25"/pixel)
	sizeAS=size
	size=size*4
	
	if size>6000:
		raise Exception('maximum size is 25 arcmin')
	
	if re.match('^[grizy]+$', band):
			pass
	else:
		raise Exception('invalid bands, supported bands: [g,r,i,z,y]')
	
	return ra,dec,size,sizeAS,band

def get_plot(ra,dec,size=30,band='g'):
	# Verify and transform inputs
	ra,dec,size,sizeAS,band=input_management(ra=ra,dec=dec,size=size,band=band)

	# Create axes and plot image
	plot=figure(width=400,height=400,title=f'PanSTARRS Image ({sizeAS}")')
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

def PanSTARRSQueryCoords(ObjRa,ObjDec,radius=3):
	radius=radius/3600

	#Had to use JSON format to get this to work, the CSV format used is not easily formattable by pandas
	url=f'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?CAT=PS1dr2OBJECTS&RA={ObjRa}&DEC={ObjDec}&SR={radius}&FORMAT=json'
	
	try:
		r=requests.get(url)
	except:
		print('[PanSTARRS: PanSTARRSQueryCoords] Error: experiencing issues with Pan-STARRS')
		return None
	
	if r.status_code!=200:
		print('[PanSTARRS: PanSTARRSQueryCoords] Error: experiencing issues with Pan-STARRS')
		return None

	data=BytesIO(r.content)
	try:
		data=pd.read_json(data)
	except:
		print('[PanSTARRS: PanSTARRSQueryCoords] Note: No data found using given coordinates')
	
	return data
	
def PanSTARRSGetPhotometryCoords(ObjRa,ObjDec,radius=3):
	data=PanSTARRSQueryCoords(ObjRa,ObjDec,radius=radius)
	
	if not isinstance(data,pd.DataFrame) or data.empty:
		print('[PanSTARRS: PanSTARRSGetPhotometryCoords] Note: No PanSTARRS data at given coordinates.')
		return None
	
	photometry=data[['ra','dec','objID','gMeanPSFMag','gMeanPSFMagErr','rMeanPSFMag','rMeanPSFMagErr','iMeanPSFMag','iMeanPSFMagErr','zMeanPSFMag','zMeanPSFMagErr','yMeanPSFMag','yMeanPSFMagErr']].copy()
	return photometry