from bokeh.plotting import figure, output_file
from bokeh.layouts import gridplot, row, column
from bokeh.models import CustomJS, Button
from bokeh import events
from bokeh.layouts import layout
import pandas as pd

'''
Changelog:
- Re-wrote everything in Bokeh
- Cleaned up imports

Known Issues:
- Detections at large distance from focus are slightly innacurate due to lack of projection support in Bokeh

To Do:
- Clean up errors
- Sort structure https://docs.python.org/3/reference/import.html#package-relative-imports
'''

# Data Queries ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def panstarrsquery(source=None,pos=None,radius=3):
	from .Surveys.PanSTARRS import PanSTARRSQueryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2012,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None

	data=PanSTARRSQueryCoords(ra,dec,radius)
	return data

def skymapperquery(source=None,pos=None,radius=3):
	from .Surveys.SkyMapper import SkyMapperQueryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2016,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=SkyMapperQueryCoords(ra,dec,radius)
	return data

def gaiaquery(source=None,pos=None,radius=3,catalogue='dr3'):
	from .Surveys.Gaia import GaiaQueryDesignation
	from .Surveys.Gaia import GaiaQueryCoords

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		data=GaiaQueryDesignation(Designation=source)
	elif pos!=None:
		ra,dec=pos[0],pos[1]
		data=GaiaQueryCoords(ra,dec)
	else:
		raise Exception('either source or pos input required')
	
	return data

def galexquery(source=None,pos=None,radius=3):
	from .Surveys.GALEX import GALEXQueryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None

	data=GALEXQueryCoords(ra,dec,radius)
	return data

def rosatquery(source=None,pos=None,radius=3):
	from .Surveys.ROSAT import ROSATQueryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[1991,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=ROSATQueryCoords(ra,dec,radius)
	return data

def sdssquery(source=None,pos=None,radius=3):
	from .Surveys.SDSS import get_data
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2017,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
		
	data=get_data(ra=ra,dec=dec,radius=radius)
	return data

def wisequery(source=None,pos=None,radius=3):
	from .Surveys.WISE import get_data
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM		

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2010,5],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
		
	data=get_data(ra=ra,dec=dec,radius=radius)
	return data

def twomassquery(source=None,pos=None,radius=3):
	from .Surveys.TWOMASS import get_data
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2010,5],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
		
	data=get_data(ra=ra,dec=dec,radius=radius)
	return data

# Imaging Queries ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getpanstarrsimage(source=None,pos=None,image_size=30,band='g',overlay=['gaia'],get_time=False):
	from astropy.time import Time

	from .Surveys.PanSTARRS import get_info
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	from .Surveys.PanSTARRS import get_plot
	from .Overlays.Overlay_Selection import overlaySelection

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		# Fetch coordinates and proper motion for object
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		# Get an image and get the time it was taken
		mjd=get_info(ra=ra,dec=dec,size=image_size,band=band)[1]
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		
		# Correct for proper motion to this image time
		pos_corrected=CorrectPM([2016,0],imageTime,ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	
	elif pos!=None:
		# No proper motion correction
		ra,dec=pos[0],pos[1]
		mjd=get_info(ra=ra,dec=dec,size=image_size,band=band)[1]
		
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
	else:
		raise Exception('either source or pos input required')
	
	# Fetch final image using coordinates corrected to the time of the original image
	plot=get_plot(ra=ra,dec=dec,size=image_size,band=band)

	if plot==None:
		if get_time==False:
			return None
		else:
			return None, None

	# Get half image size (in deg, used for detection size scaling)
	border=image_size/7200

	if source!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border,pmra,pmdec)
	elif pos!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border)

	if plot!=None and detections_made!=False:
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

	if source!=None:
		output_file(f'{source}_image.html')

	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_image.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_image.html")

	if get_time==True:
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		return plot, imageTime
	else:
		return plot

def getskymapperimage(source=None,pos=None,image_size=30,band='g',overlay=['gaia'],get_time=False):
	from astropy.time import Time
	
	from .Surveys.SkyMapper import get_info
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	from .Surveys.SkyMapper import get_plot
	from .Overlays.Overlay_Selection import overlaySelection

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		# Fetch coordinates and proper motion for object
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		# Get an image and get the time it was taken
		mjd=get_info(ra=ra,dec=dec,size=image_size,band=band)[1]
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		
		# Correct for proper motion to this image time
		pos_corrected=CorrectPM([2016,0],imageTime,ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	
	elif pos!=None:
		# No proper motion correction
		ra,dec=pos[0],pos[1]
		mjd=get_info(ra=ra,dec=dec,size=image_size,band=band)[1]
		
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
	else:
		raise Exception('either source or pos input required')
	
	# Fetch final image using coordinates corrected to the time of the original image
	plot=get_plot(ra=ra,dec=dec,size=image_size,band=band)

	if plot==None:
		if get_time==False:
			return None
		else:
			return None, None

	# Get half image size (in deg, used for detection size scaling)
	border=image_size/7200

	if source!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border,pmra,pmdec)
	elif pos!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border)
	
	if plot!=None and detections_made!=False:
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

	if source!=None:
		output_file(f'{source}_image.html')

	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_image.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_image.html")

	if get_time==True:
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		return plot, imageTime
	else:
		return plot

def getdssimage(source=None,pos=None,image_size=30,overlay=['gaia'],get_time=False):
	from astropy.time import Time
	
	from .Surveys.DSS import get_info
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	from .Surveys.DSS import get_plot
	from .Overlays.Overlay_Selection import overlaySelection

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		# Fetch coordinates and proper motion for object
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		# Get an image and get the time it was taken
		mjd=get_info(ra=ra,dec=dec,size=image_size)[1]
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
			
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		
		# Correct for proper motion to this image time
		pos_corrected=CorrectPM([2016,0],imageTime,ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	
	elif pos!=None:
		# No proper motion correction
		ra,dec=pos[0],pos[1]
		mjd=get_info(ra=ra,dec=dec,size=image_size)[1]
			
		if mjd==None:
			if get_time==False:
				return None
			else:
				return None, None
	else:
		raise Exception('either source or pos input required')
	
	# Fetch final image using coordinates corrected to the time of the original image
	plot=get_plot(ra=ra,dec=dec,size=image_size)

	if plot==None:
		if get_time==False:
			return None
		else:
			return None, None

	# Get half image size (in deg, used for detection size scaling)
	border=image_size/7200

	if source!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border,pmra,pmdec)
	elif pos!=None:
		plot,detections_made=overlaySelection(plot,ra,dec,overlay,mjd,image_size,border)

	plot.legend.click_policy="hide"	

	if plot!=None and detections_made!=None:
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

	if source!=None:
		output_file(f'{source}_image.html')
	
	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_image.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_image.html")

	if get_time==True:
		imageTime=Time(mjd,format='mjd').to_datetime()
		imageTime=[imageTime.year,imageTime.month]
		return plot, imageTime
	else:
		return plot

# Photometry Queries  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getpanstarrsphot(radius=3,source=None,pos=None):
	from .Surveys.PanSTARRS import PanSTARRSGetPhotometryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2012,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
		
	data=PanSTARRSGetPhotometryCoords(ra,dec,radius)
	return data

def getskymapperphot(radius=3,source=None,pos=None):
	from .Surveys.SkyMapper import SkyMapperGetPhotometryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2016,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None

	data=SkyMapperGetPhotometryCoords(ra,dec,radius)
	return data

def getgaiaphot(radius=3,source=None,pos=None):
	from .Surveys.Gaia import GaiaGetPhotometryCoords
	from .Surveys.Gaia import GaiaGetPhotometryDesignation
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		data=GaiaGetPhotometryDesignation(Designation=source)
	elif pos!=None:
		ra,dec=pos[0],pos[1]
		data=GaiaGetPhotometryCoords(ra,dec,radius)
	else:
		raise Exception('either source or pos input required')
		return None

	return data

def getgalexphot(radius=3,source=None,pos=None):
	from .Surveys.GALEX import GALEXGetPhotometryCoords
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None

	data=GALEXGetPhotometryCoords(ra,dec,radius)
	
	return data

def getsdssphot(radius=3,source=None,pos=None):
	from .Surveys.SDSS import get_photometry
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=get_photometry(ra=ra,dec=dec,radius=radius)
	return data

def getwisephot(radius=3,source=None,pos=None):
	from .Surveys.WISE import get_photometry
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=get_photometry(ra=ra,dec=dec,radius=radius)
	return data

def gettwomassphot(radius=3,source=None,pos=None):
	from .Surveys.TWOMASS import get_photometry
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM	
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=get_photometry(ra=ra,dec=dec,radius=radius)
	return data

# Bulk Photometry Query

def getbulkphot(radius=3,source=None,pos=None):
	from .Surveys.Gaia import GaiaGetPhotometryCoords
	from .Surveys.GALEX import GALEXGetPhotometryCoords
	from .Surveys.ROSAT import ROSATGetPhotometryCoords
	from .Surveys.PanSTARRS import PanSTARRSGetPhotometryCoords
	from .Surveys.SkyMapper import SkyMapperGetPhotometryCoords
	from .Surveys.SDSS import get_photometry as get_phot_sdss
	from .Surveys.WISE import get_photometry as get_phot_wise
	from .Surveys.TWOMASS import get_photometry as get_phot_twomass
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:		
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2012,0],ra,dec,pmra,pmdec)
		ra_panstarrs,dec_panstarrs=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=CorrectPM([2016,0],[2016,0],ra,dec,pmra,pmdec)
		ra_skymapper,dec_skymapper=pos_corrected[0],pos_corrected[1]
		
		ra_gaia,dec_gaia=ra,dec
		
		pos_corrected=CorrectPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra_galex,dec_galex=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=CorrectPM([2016,0],[1991,0],ra,dec,pmra,pmdec)
		ra_rosat,dec_rosat=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=CorrectPM([2016,0],[2017,0],ra,dec,pmra,pmdec)
		ra_sdss,dec_sdss=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=CorrectPM([2016,0],[2010,5],ra,dec,pmra,pmdec)
		ra_wise,dec_wise=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=CorrectPM([2016,0],[1999,0],ra,dec,pmra,pmdec)
		ra_twomass,dec_twomass=pos_corrected[0],pos_corrected[1]
		
	elif pos!=None:
		ra_panstarrs,dec_panstarrs=pos[0],pos[1]
		ra_skymapper,dec_skymapper=pos[0],pos[1]
		ra_gaia,dec_gaia=pos[0],pos[1]
		ra_galex,dec_galex=pos[0],pos[1]
		ra_rosat,dec_rosat=pos[0],pos[1]
		ra_sdss,dec_sdss=pos[0],pos[1]
		ra_wise,dec_wise=pos[0],pos[1]
		ra_twomass,dec_twomass=pos[0],pos[1]
		
	else:
		raise Exception('either source or pos input required')
		return None

	photometry={'gaia':None,'galex':None,'rosat':None,'panstarrs':None,'skymapper':None,'sdss':None,'wise':None,'twomass':None}
	
	try:
		data=GaiaGetPhotometryCoords(ra_gaia,dec_gaia,radius)
		photometry['gaia']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No Gaia photometry found using given coordinates and search radius')
	try:
		data=GALEXGetPhotometryCoords(ra_galex,dec_galex,radius)
		photometry['galex']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No GALEX photometry found using given coordinates and search radius')
	try:
		data=ROSATGetPhotometryCoords(ra_rosat,dec_rosat,radius)
		photometry['rosat']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No ROSAT photometry found using given coordinates and search radius')
	try:
		data=PanSTARRSGetPhotometryCoords(ra_panstarrs,dec_panstarrs,radius)
		photometry['panstarrs']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No Pan-STARRS photometry found using given coordinates and search radius')
	try:
		data=SkyMapperGetPhotometryCoords(ra_skymapper,dec_skymapper,radius)
		photometry['skymapper']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No SkyMapper photometry found using given coordinates and search radius')
	try:
		data=get_phot_sdss(ra=ra_sdss,dec=dec_sdss,radius=radius)
		photometry['sdss']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No SDSS photometry found using given coordinates and search radius')
	try:
		data=get_phot_wise(ra=ra_wise,dec=dec_wise,radius=radius)
		photometry['wise']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No WISE photometry found using given coordinates and search radius')
	try:
		data=get_phot_twomass(ra=ra_twomass,dec=dec_twomass,radius=radius)
		photometry['twomass']=data
	except:
		print('[Photometry: GetPhotometryCoords] Note: No 2MASS photometry found using given coordinates and search radius')
	
	return photometry

# Timeseries Queries -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def ztfquery(source=None,pos=None,radius=3):
	from .Surveys.ZTF import getData as getZTFData	
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2019,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	data=getZTFData(ra,dec,radius)
	return data

# Timeseries Plotting ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getztflc(source=None,pos=None,radius=3,return_raw=False):
	from .Surveys.ZTF import getLightCurve as getZTFLightCurve
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2019,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
		
	plot=getZTFLightCurve(ra,dec,radius,return_raw)
	
	if source!=None:
		output_file(f'{source}_lightcurve.html')
	
	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_lightcurve.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_lightcurve.html")

	return plot

# SED Plotting -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getsed(source=None,pos=None,radius=3):
	from .Figures.SED import get_plot
	
	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')

	if source!=None:
		plot=get_plot(source=source,radius=radius)
	elif pos!=None:
		plot=get_plot(pos=pos,radius=radius)
	else:
		raise Exception('either source or pos input required')

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
	else:
		return None

	if source!=None:
		output_file(f'{source}_sed.html')

	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_sed.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_sed.html")
		
	return plot

# Spectra Queries --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getsdssspectrum(source=None,pos=None,radius=3):
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM
	from .Surveys.SDSS import get_plot

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=CorrectPM([2016,0],[2017,0],ra,dec,pmra,pmdec)
		ra,dec=pos_corrected[0],pos_corrected[1]
	elif pos!=None:
		ra,dec=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')
		return None
	
	plot=get_plot(ra=ra,dec=dec,radius=radius)

	if source!=None:
		output_file(f'{source}_spectrum.html')

	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_spectrum.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_spectrum.html")
		
	return plot

# HR diagram --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def gethrd(source=None,sources=None):
	from .Figures.HRD import get_plot
	if source==None and sources==None:
		raise Exception('source/sources input required')
		return None
	
	if source!=None:
		plot=get_plot(source=source)
	elif sources!=None:
		plot=get_plot(sources=sources)
	else:
		raise Exception('source/sources input required.')
	
	if plot!=None:
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

	if source!=None:
		output_file(f'{source}_hrd.html')
	elif sources==None:
		sources_str=''
		for i in range(0,len(sources)-1):
			sources_str.append(str(sources[i]))+','
		sources_str.append(source[len(source)-1])
		output_file(f'{sources}_hrd.html')

	return plot

# Timeseries analysis -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getztfanalysis(source=None,pos=None):
	from .Timeseries.ztfanalysis import getanalysis
	
	if pos!=None:
		data=ztfquery(pos=pos)
	elif source!=None:
		data=ztfquery(source=source)
	else:
		raise Exception('either source or pos input required')
		return None

	empty_count=0
	for item in data:
		if not isinstance(item,pd.DataFrame):
			empty_count+=1
			
	if empty_count!=3:
		data=pd.concat(data)
	else:
		print('no ZTF data available for given fields')
		return None

	getanalysis(data)

def getps(source=None,pos=None):
	from .Timeseries.ztfanalysis import getpowerspectrum
	
	if pos!=None:
		data=ztfquery(pos=pos)
	elif source!=None:
		data=ztfquery(source=source)
	else:
		raise Exception('either source or pos input required')

	empty_count=0
	for item in data:
		if not isinstance(item,pd.DataFrame):
			empty_count+=1
	
	if empty_count!=3:
		data=pd.concat(data)
	else:
		print('no ZTF data available for given fields')
		return None	

	plot=getpowerspectrum(data)

	if source!=None:
		output_file(f'{source}_powspec.html')

	elif pos!=None:
		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_powspec.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_powspec.html")

	return plot

# Miscellaneous Tools -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def correctPM(input,target,ra,dec,pmra,pmdec):
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	data=CorrectPM(input,target,ra,dec,pmra,pmdec)
	return data

def getgaiacoords(source,catalogue='dr3'):
	from .Surveys.Gaia import GaiaGetCoords
	
	data=GaiaGetCoords(source,catalogue)
	return data

def getgaiasource(pos,radius=3,catalogue='dr3'):
	from .Surveys.Gaia import GaiaGetDesignation
	
	ra,dec=pos[0],pos[1]
	data=GaiaGetDesignation(ra,dec,radius,catalogue)
	return data

def getsources(file_name):
	from .Miscellaneous.ReadFits import get_source_list
	
	sources=get_source_list(file_name)
	return sources

def getpositions(file_name):
	from .Miscellaneous.ReadFits import get_pos_list
	pos_list=get_pos_list(file_name)
	return pos_list

# Datapage generation -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getdatapage(source=None,pos=None,prefs=None):
	from .Miscellaneous.ProperMotionCorrection import PMCorrection as CorrectPM

	if isinstance(prefs,dict):
		keys=prefs.keys()
		if 'band' in keys:
			band=prefs['band']
		else:
			band='g'
		if 'overlay' in keys:
			overlay=prefs['overlay']
		else:
			overlay=['gaia']
		if 'image_size' in keys:
			image_size=prefs['image_size']
		else:
			image_size=20
		if 'lc_radius' in keys:
			lc_radius=prefs['lc_radius']
		else:
			lc_radius=3
		if 'sed_radius' in keys:
			sed_radius=prefs['sed_radius']
		else:
			sed_radius=3
		if 'spec_radius' in keys:
			spec_radius=prefs['spec_radius']
		else:
			spec_radius=3
		if 'simbad_radius' in keys:
			simbad_radius=prefs['simbad_radius']
		else:
			simbad_radius=3
		if 'vizier_radius' in keys:
			vizier_radius_override=prefs['vizier_radius']
		else:
			vizier_radius_override=None
		if 'grid_size' in keys:
			grid_size_override=prefs['grid_size']
		else:
			grid_size_override=300

	else:
		band='g'
		overlay=['gaia']
		image_size=20
		lc_radius=3
		sed_radius=3
		spec_radius=3
		simbad_radius=3
		vizier_radius_override=None
		grid_size_override=300

	if source!=None and pos!=None:
		raise Exception('simulatenous source and pos input detected')
	
	if source!=None:
		gaia_data=gaiaquery(source=source)
		ra,dec,pmra,pmdec=gaia_data['ra'].values[0],gaia_data['dec'].values[0],gaia_data['pmra'].values[0],gaia_data['pmdec'].values[0]
		
		pos_corrected=correctPM([2016,0],[2012,0],ra,dec,pmra,pmdec)
		ra_panstarrs,dec_panstarrs=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=correctPM([2016,0],[2016,0],ra,dec,pmra,pmdec)
		ra_skymapper,dec_skymapper=pos_corrected[0],pos_corrected[1]
		
		ra_gaia,dec_gaia=ra,dec
		
		pos_corrected=correctPM([2016,0],[2007,0],ra,dec,pmra,pmdec)
		ra_galex,dec_galex=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=correctPM([2016,0],[1991,0],ra,dec,pmra,pmdec)
		ra_rosat,dec_rosat=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=correctPM([2016,0],[2019,0],ra,dec,pmra,pmdec)
		ra_ztf,dec_ztf=pos_corrected[0],pos_corrected[1]
		
		pos_corrected=correctPM([2016,0],[2017,0],ra,dec,pmra,pmdec)
		ra_sdss,dec_sdss=pos_corrected[0],pos_corrected[1]
		
	elif pos!=None:
		ra,dec=pos[0],pos[1]
		ra_panstarrs,dec_panstarrs=pos[0],pos[1]
		ra_skymapper,dec_skymapper=pos[0],pos[1]
		ra_gaia,dec_gaia=pos[0],pos[1]
		ra_galex,dec_galex=pos[0],pos[1]
		ra_rosat,dec_rosat=pos[0],pos[1]
		ra_ztf,dec_ztf=pos[0],pos[1]
		ra_sdss,dec_sdss=pos[0],pos[1]
	else:
		raise Exception('either source or pos input required')

	# Need to add in 'prefs' parameters
	if source!=None:
		lightcurve_axis_g,lightcurve_axis_r,lightcurve_axis_i=getztflc(source=source,return_raw=True,radius=lc_radius)
		hrd_axis=gethrd(source=source)
		sed_axis=getsed(source=source,radius=sed_radius)
		spectrum_axis=getsdssspectrum(source=source,radius=spec_radius)
		image_axis,image_time=getpanstarrsimage(source=source,get_time=True,band=band,overlay=overlay,image_size=image_size)
		if image_axis==None:
			image_axis,image_time=getskymapperimage(source=source,get_time=True,band=band,overlay=overlay,image_size=image_size)
			if image_axis==None:
				image_axis,image_time=getdssimage(source=source,get_time=True,overlay=overlay,image_size=image_size)
		ps_axis=getps(source=source)

		output_file(f"{source}_datapage.html")	
	
	elif pos!=None:
		lightcurve_axis_g,lightcurve_axis_r,lightcurve_axis_i=getztflc(pos=pos,return_raw=True,radius=lc_radius)
		hrd_axis=None
		sed_axis=getsed(pos=pos,radius=sed_radius)
		spectrum_axis=getsdssspectrum(pos=pos,radius=spec_radius)
		image_axis,image_time=getpanstarrsimage(pos=pos,get_time=True,band=band,overlay=overlay,image_size=image_size)
		if image_axis==None:
			image_axis,image_time=getskymapperimage(pos=pos,get_time=True,band=band,overlay=overlay,image_size=image_size)
			if image_axis==None:
				image_axis,image_time=getdssimage(pos=pos,get_time=True,overlay=overlay,image_size=image_size)
		ps_axis=getps(pos=pos)

		if pos[1]>=0:
			output_file(f"{pos[0]}+{pos[1]}_datapage.html")	
		else:
			output_file(f"{pos[0]}{pos[1]}_datapage.html")	

	button_width=100
	button_height=25

	# SIMBAD button
	simbad_button = Button(label="SIMBAD",button_type='primary',height=button_height,width=button_width)	

	if pos!=None:
		simbad_url=f'https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={ra}+{dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={simbad_radius}&Radius.unit=arcsec&submit=submit+query&CoordList='
	elif source!=None:
		simbad_url=f'http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Gaia DR3 {source}&NbIdent=1&Radius={simbad_radius}&Radius.unit=arcsec&submit=submit+id'
	
	simbad_button_js = CustomJS(args=dict(url=simbad_url),code='''
		window.open(url)
	''')
	simbad_button.js_on_event('button_click',simbad_button_js)

	# Vizier button
	# Scale vizier search radius and coordinates with proper motion if available
	if source!=None:
		if isinstance(image_time,list):
			gaia_time=[2016,0]
			vizier_ra,vizier_dec,vizier_radius=CorrectPM(gaia_time,image_time,ra_gaia,dec_gaia,pmra,pmdec,radius=image_size)
		else:
			vizier_ra,vizier_dec,vizier_radius=ra,dec,image_size
	elif pos!=None:
		vizier_ra,vizier_dec,vizier_radius=pos[0],pos[1],image_size

	# Allows search radius to be overriden by prefs
	if vizier_radius_override!=None:
		vizier_radius=vizier_radius_override

	vizier_button = Button(label="Vizier",button_type='primary',height=button_height,width=button_width)	
	
	if vizier_dec>=0:
		vizier_url=f'https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={vizier_ra}+{vizier_dec}&-c.rs={vizier_radius}&-out.add=_r&-sort=_r&-out.max=$4'
	else:
		vizier_url=f'https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={vizier_ra}{vizier_dec}&-c.rs={vizier_radius}&-out.add=_r&-sort=_r&-out.max=$4'
	
	vizier_button_js = CustomJS(args=dict(url=vizier_url),code='''
		window.open(url)
	''')
	vizier_button.js_on_event('button_click',vizier_button_js)

	# Set up grid and button dimensions
	grid_size=grid_size_override
	
	# Empirical relation that roughly works
	left_offset=round((0.25+((300-grid_size)/1000))*grid_size)
	
	# Margin = [Top,Right,Bottom,Left]
	simbad_button.margin=[round(0.1*(grid_size-button_height)),0,0,left_offset]
	vizier_button.margin=[round(0.05*(grid_size-button_height)),0,0,left_offset]

	simbad_button.width,simbad_button.height=round(0.6*grid_size),round(0.15*grid_size)
	vizier_button.width,vizier_button.height=round(0.6*grid_size),round(0.15*grid_size)

	# layout doesn't allow None
	if image_axis==None:
		image_axis=figure()
	if lightcurve_axis_g==None:
		lightcurve_axis_g=figure()
	if lightcurve_axis_r==None:
		lightcurve_axis_r=figure()
	if lightcurve_axis_i==None:
		lightcurve_axis_i=figure()
	if sed_axis==None:
		sed_axis=figure()
	if hrd_axis==None:
		hrd_axis=figure()
	if spectrum_axis==None:
		spectrum_axis=figure()
	if ps_axis==None:
		ps_axis=figure()

	buttons=column(simbad_button,vizier_button,align='center')

	image_axis.width,image_axis.height=2*grid_size,2*grid_size
	buttons.width,buttons.height=grid_size,grid_size
	lightcurve_axis_g.width,lightcurve_axis_g.height=grid_size,grid_size
	
	hrd_axis.width,hrd_axis.height=2*grid_size,2*grid_size
	lightcurve_axis_r.width,lightcurve_axis_r.height=grid_size,grid_size
	lightcurve_axis_i.width,lightcurve_axis_i.height=grid_size,grid_size

	sed_axis.width,sed_axis.height=2*grid_size,grid_size
	spectrum_axis.width,spectrum_axis.height=2*grid_size,grid_size
	ps_axis.width,ps_axis.height=2*grid_size,grid_size

	grid=layout(row(
		column([image_axis,row(buttons,lightcurve_axis_g)]),
		column([hrd_axis,row(lightcurve_axis_r,lightcurve_axis_i)]),
		column([sed_axis,spectrum_axis,ps_axis])
	))
	
	'''
	def exec_timeseries(event):
		sys_path=os.getcwd()
		timeseries_path=os.path.join(sys_path,'timeseries.py')
		print(timeseries_path)
		if source!=None:
			subprocess.run(['python',timeseries_path,'source',str(source)])
		elif pos!=None:
			subprocess.run(['python',timeseries_path,'pos',str(pos[0]),str(pos[1])])
	
	timeseries_button.on_event('button_click',exec_timeseries)
	'''

	return grid