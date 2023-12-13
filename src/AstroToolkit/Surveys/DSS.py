from astropy.io import fits
from astropy.visualization import PercentileInterval, AsinhStretch
from io import BytesIO
import requests
import numpy as np
from astropy.wcs import WCS
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.time import Time
import astropy
from bokeh.plotting import figure
from bokeh.models import Range1d


#Ignore warning about deprecated header
warnings.simplefilter('ignore',category=AstropyWarning)

# DSS Query
def get_image(ra,dec,size=30):
	service='http://archive.stsci.edu/cgi-bin/dss_search'
	url=f'{service}?ra={ra}&d={dec}&v=3&e=J2000&f=fits&h={size}&w={size}'
	
	try:
		r=requests.get(url,timeout=15)
		if r.status_code!=200:
			print('[DSS: get_image_data] Error: experiencing issues with DSS.')
			return None, None
	except:
		print('[DSS: get_image_data] Error: experiencing issues with DSS.')
		return None, None
	
	try:
		fh=fits.open(BytesIO(r.content))
		fitsImage=fh[0].data
		image_header=fh[0].header
		fitsImage[np.isnan(fitsImage)]=0.0
		transform=AsinhStretch()+PercentileInterval(95)
		image=transform(fitsImage)
		
		if not isinstance(image,np.ndarray):
			return None, None

		return image,image_header
	except:
		print('[DSS: get_image_data] Error: no DSS images found using given fields')
		return None, None

# Filters inputs and performs any necessary arithmetic
def input_management(ra,dec,size=30):
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
	
	if not isinstance(size,int):
		try:
			size=int(size)
		except:
			raise Exception('size must take an integer value')
	
	sizeAS=size
	size=size/60
	
	if size>60:
		raise Exception('maximum size allowed by DSS is 60 arcminutes.')
	
	return ra, dec, size, sizeAS

# Returns only the wcs projection for a retrieved image
def get_info(ra,dec,size=30):
	ra,dec,size,sizeAS=input_management(ra=ra,dec=dec,size=size)

	image_data=get_image(ra=ra,dec=dec,size=size)[1]
	
	if not isinstance(image_data,astropy.io.fits.header.Header):
		return None, None
	
	wcs=WCS(image_data)
	
	time=image_data['DATE-OBS']
	mins=int(time[14:16])
	hours=int(time[11:13])
	if mins>=60:
		mins=mins-60
		hours+=1
		
	mins=str(mins).zfill(2)
	hours=str(hours).zfill(2)
	
	time=time[0:10]+'T'+str(hours)+':'+str(mins)+':'+time[17:20]
	
	imageTime=Time(time,format='fits')
	
	mjd=imageTime.mjd
	
	return wcs, mjd

# Sets up axis and casts image onto them
def get_plot(ra,dec,size=30):
	ra,dec,size,sizeAS=input_management(ra=ra,dec=dec,size=size)

	# Create axes and plot image
	plot=figure(width=400,height=400,title=f'DSS Image ({sizeAS}")')
	# Fixes scaling issue
	plot.min_border=75

	plot.xgrid.grid_line_color=None
	plot.ygrid.grid_line_color=None	

	image,image_data=get_image(ra=ra,dec=dec,size=size)
	
	if not isinstance(image,np.ndarray) or not isinstance(image_data,astropy.io.fits.header.Header):
		return None

	# Get x and y limit (in pixels) from image header
	xlim=image_data['NAXIS1']
	ylim=image_data['NAXIS2']

	# DSS images have one less pixel in one axis sometimes, so set both axes limits to that of the smaller one
	lower_limit=min([xlim,ylim])
	if lower_limit==xlim:
		ylim=xlim
	elif lower_limit==ylim:
		xlim=ylim
	
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