from tkinter import W
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
import warnings
from astropy.utils.exceptions import AstropyWarning
import pandas as pd
from astropy.time import Time

from ..Data.data import get_survey_times
survey_times=get_survey_times()

# ignores a warning about a deprecated header
warnings.simplefilter('ignore', category=AstropyWarning)

'''
Fetches image data dictionary
'''
def get_image_dict(survey,source=None,pos=None,location=None,size=30,band='g'):
	def failed_search():
		return {'type':'image','survey':survey,'source':source,'pos':pos,'data':None}
	
	ra,dec=location[0],location[1]

	if survey=='panstarrs':
		# convert size to pixel size
		url_size=size*4
		service='https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
		url=f'{service}?ra={ra}&dec={dec}&band={band}'
	
		# get table of returned image filenames
		try:
			table=Table.read(url, format='ascii')
		except:
			print(f'Note: Experiencing issues with {survey}.')
			return failed_search()
		
		if not len(table)>0:
			print(f'Note: No data returned from {survey} imagequery.')
			return failed_search()
		
		url=(f'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={ra}&dec={dec}&size={url_size}&format=fits')
	
		filterList=['yzirg'.find(x) for x in table['filter']]
		table=table[np.argsort(filterList)]
	
		urlbase=url+'&red='
		url=[]
		for filename in table['filename']:
			url.append(urlbase+filename)
		
		# return URL for first image from the above table
		urlMain=url[0]
		
	elif survey=='skymapper':
		# convert size to pixel size
		url_size=size/3600

		service='https://api.skymapper.nci.org.au/public/siap/dr2/query'
		url=f'{service}?POS={ra},{dec}&SIZE={url_size},{url_size}&BAND={band}&FORMAT=image/fits&VERB=3&INTERSECT=covers&RESPONSEFORMAT=CSV'
		
		# get table of available images
		try:
			table=pd.read_csv(url)
			table=Table.from_pandas(table)
		except:
			print(f'Note: Experiencing issues with {survey}.')
			return failed_search()
	
		# get url of first image from the above table
		if len(table)>0:
			urlMain=table['get_image'][0]
		else:
			print(f'Note: No data returned from {survey} imagequery.')
			return failed_search()
	
	elif survey=='dss':
		# convert size to pixel size
		url_size=size/60
		
		service='http://archive.stsci.edu/cgi-bin/dss_search'
		urlMain=f'{service}?ra={ra}&d={dec}&v=3&e=J2000&f=fits&h={url_size}&w={url_size}'

	# send request using the URL returned for the selected survey
	r=requests.get(urlMain,timeout=15)
	if r.status_code!=200:
		print(f'Note: Experiencing issues with {survey}.')

	try:
		fh=fits.open(BytesIO(r.content))
	except:
		print(f'Note: No data returned from {survey} imagequery.')
		return failed_search()
	
	# read image fits file, apply a contrast filter
	fitsImage=fh[0].data
	image_header=fh[0].header
	fitsImage[np.isnan(fitsImage)]=0.0
	transform=AsinhStretch()+PercentileInterval(95)
	image=transform(fitsImage)
	
	# get imageTime [year,month] of image from its header
	if survey=='panstarrs' or survey=='skymapper':
		mjd=image_header['MJD-OBS']
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
	
	# dss is not as simple, and sometimes returns erroneous values (i.e. minutes>60), so need to fix this
	elif survey=='dss':
		time=image_header['DATE-OBS']
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
		
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]

	# get world coordinates system used to turn pixel coordinates into world coordinates
	wcs=WCS(image_header)
	
	# setup image dictionary
	image_dict={'type':'image','survey':survey,'source':source,'pos':pos,'data':{'image_data':image,'image_header':image_header,'location':[ra,dec],'size':size,'image_time':imageTime,'wcs':wcs,'overlay':None}}
		
	return image_dict

'''
returns image dictionaries after correcting for proper motion (if source is used as input)
'''
def image_correction(survey,source=None,pos=None,size=30,band='g',overlay=None):
	from ..Tools import dataquery
	from ..Misc.ProperMotionCorrection import PMCorrection
	from ..Data.overlays import overlay_selection

	if source!=None:
		# get Gaia data, if None is returned this will fail, triggering except --> returns None

		gaia_data=dataquery(survey='gaia',source=source)['data']
		if gaia_data!=None:
			ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]
		else:
			raise Exception('Gaia source not found.')

		# get an initial image, which is just used to see when the image was taken
		image_dict=get_image_dict(survey=survey,pos=pos,source=source,location=[ra,dec],size=size,band=band)

		# get initial image time
		if image_dict['data']!=None:
			image_time=image_dict['data']['image_time']
		else:
			return image_dict

		# correct source's coordinates for proper motion back to this image time
		pos_corrected=PMCorrection(survey_times['gaia'],image_time,ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
		
		# get another image dict at these new coordinates, and apply overlay if an image is returned
		image_dict=get_image_dict(survey=survey,pos=pos,source=source,location=[ra,dec],size=size,band=band)
		if image_dict['data']!=None:
			image_dict=overlay_selection(image_dict,overlay)
		else:
			return image_dict

	# no proper motion correction, just take an image and apply an overlay
	elif pos!=None:
		ra,dec=pos[0],pos[1]
		image_dict=get_image_dict(survey=survey,source=source,pos=pos,location=[ra,dec],size=size,band=band)
		
		if image_dict['data']!=None:
			image_dict=overlay_selection(image_dict=image_dict,overlay=overlay)
		else:
			return image_dict
		
	return image_dict