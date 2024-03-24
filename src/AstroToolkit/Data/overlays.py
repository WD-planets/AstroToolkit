import pandas as pd
import math

from ..Misc.ProperMotionCorrection import PMCorrection
from ..Misc.ProperMotionCorrection import get_adapted_radius
from ..Tools import dataquery
from ..Tools import lightcurvequery
from ..Data.data import get_survey_times

survey_times=get_survey_times()

'''
get overlay data for a given image and survey. overlay data is in format [{overlay_entry},{overlay_entry},...] with each overlay entry giving the information of a single overlay detection
'''
def get_overlay(image_dict,survey,ra_name,dec_name,mag_name,overlay_type='detection',marker_type='circle'):
	# 'detections' here are any non-timeseries surveys (i.e. just point to show that data exists)
	if overlay_type=='detection':
		# radius used for cross-matching with gaia detections to correct for proper motion in other surveys
		piggyback_radius=5
	
		gaia_time=[2016,0]
		source=image_dict['source']
		location=image_dict['data']['location']
		size=image_dict['data']['size']
		half_image=image_dict['data']['size']/7200
		image_time=image_dict['data']['image_time']

		# left here as the size of the image isn't always the right radius to use in a query (some use circular, some use square), etc.
		radius=size

		# get gaia data for the source object, and use its proper motion to scale the radius used.
		if source!=None:
			gaia_data=dataquery(survey='gaia',source=source)['data']
			if gaia_data!=None:
				pmra,pmdec=gaia_data['pmra'][0],gaia_data['pmdec'][0]
				radius=get_adapted_radius(input=gaia_time,target=image_time,pmra=pmra,pmdec=pmdec,radius=radius)
			else:
				return None

		# get overlay survey data
		nearby_non_gaia=dataquery(survey=survey,pos=location,radius=radius)['data']
		if nearby_non_gaia==None:
			return None
		else:
			non_gaia_count=len(nearby_non_gaia[list(nearby_non_gaia.keys())[0]])
	
		# get nearby gaia data (now using the radius calculated above)
		nearby_gaia=dataquery(survey='gaia',pos=location,radius=radius)['data']
		if nearby_gaia!=None:
			gaia_count=len(nearby_gaia[list(nearby_gaia.keys())[0]])
				
		skip_mag=False
		nearby_survey_ra=nearby_non_gaia[ra_name]
		nearby_survey_dec=nearby_non_gaia[dec_name]
		if mag_name!=None:
			nearby_survey_mag=nearby_non_gaia[mag_name]
			# check if any magnitudes exist in the returned non-gaia data	
			mag_check=False
			for mag in nearby_survey_mag:
				if not math.isnan(mag):
					mag_check=True
		
			if mag_check!=True:
				return None
		else:
			skip_mag=True

		# check if any gaia data was returned
		if nearby_gaia!=None:
			nearby_gaia_ra=nearby_gaia['ra']
			nearby_gaia_dec=nearby_gaia['dec']
			nearby_gaia_pmra=nearby_gaia['pmra']
			nearby_gaia_pmdec=nearby_gaia['pmdec']

			# corrects gaia coordinates for proper motion back to the time of the non-gaia survey
			for i in range(0,gaia_count):
				if not math.isnan(nearby_gaia_pmra[i]) and not math.isnan(nearby_gaia_pmdec[i]):
					pos_corrected=PMCorrection(input=survey_times['gaia'],target=survey_times[survey],ra=nearby_gaia_ra[i],dec=nearby_gaia_dec[i],pmra=nearby_gaia_pmra[i],pmdec=nearby_gaia_pmdec[i])
					nearby_gaia_ra[i],nearby_gaia_dec[i]=pos_corrected[0],pos_corrected[1]

			# checks in a radius around these corrected gaia detections. if a detection from the non-gaia survey is found within the piggyback_radius, applies the reverse of the proper motion correction that was applied to the gaia correction to this non-gaia detection
			corrected_indices=[]
			for i in range(0,non_gaia_count):
				for j in range(0,gaia_count):
					delta=math.sqrt((nearby_survey_ra[i]-nearby_gaia_ra[j])**2+(nearby_survey_dec[i]-nearby_gaia_dec[j])**2)*3600
					if delta<piggyback_radius:
						if not math.isnan(nearby_gaia_pmra[j]) and not math.isnan(nearby_gaia_pmdec[j]):
							pos_corrected=PMCorrection(input=survey_times[survey],target=image_time,ra=nearby_survey_ra[i],dec=nearby_survey_dec[i],pmra=nearby_gaia_pmra[j],pmdec=nearby_gaia_pmdec[j])
							nearby_survey_ra[i],nearby_survey_dec[i]=pos_corrected[0],pos_corrected[1]
							corrected_indices.append(i)
		
			# scales detection radius based on the given magnitude
			overlay_data=[]
			for i in range(0,non_gaia_count):
				if skip_mag!=True:
					radiusMultiplier=(nearby_survey_mag[i]/20.7)
					radius=(half_image/50+(half_image/75)**radiusMultiplier)*0.75
				else:
					radius=None
	
				if i in corrected_indices:
					corrected=True
				else:
					corrected=False

				# set ups an overlay entry, and appends it to the list of entries for this survey
				overlay_entry={'survey':survey,'position':[nearby_survey_ra[i],nearby_survey_dec[i]],'radius':radius,'corrected':corrected,'mag':mag_name,'marker':marker_type}
				overlay_data.append(overlay_entry)
		
		# if no gaia data was returned, cannot do any 'piggybacking' with other surveys, and so skips this and just marks them as uncorrected by default
		else:
			overlay_data=[]
			for i in range(0,non_gaia_count):
				if skip_mag!=True:
					radiusMultiplier=(nearby_survey_mag[i]/20.7)
					radius=(half_image/50+(half_image/75)**radiusMultiplier)*0.75
				else:
					radius=None
				
				# set ups an overlay entry, and appends it to the list of entries for this survey
				overlay_entry={'survey':survey,'position':[nearby_survey_ra[i],nearby_survey_dec[i]],'radius':radius,'corrected':False,'mag':mag_name,'marker':marker_type}
				overlay_data.append(overlay_entry)
	
	# 'tracer' here refers to timeseries surveys (e.g. ztf), which can 'trace' an object through time
	elif overlay_type=='tracer':
		source=image_dict['metadata']['source']
		location=image_dict['metadata']['location']
		size=image_dict['metadata']['size']		
			
		if survey=='ztf':
			# set ztf radius to a circle with radius = diagonal of square image
			radius=math.sqrt(2*(size/2)**2)
		elif survey=='atlas' or survey=='gaia_lc' or survey=='asassn' or survey=='crts':
			radius=size

		# get lightcurve data in format [band1,band2,etc]
		if survey!='gaia_lc':
			if source==None:
				data_arr=lightcurvequery(survey=survey,pos=location,radius=radius)
			else:
				data_arr=lightcurvequery(survey=survey,source=source,radius=radius)
		else:
			if source==None:
				data_arr=lightcurvequery(survey='gaia',pos=location,radius=radius)
			else:
				data_arr=lightcurvequery(survey='gaia',source=source,radius=radius)

		data_exists=False		
		for i,val in enumerate(data_arr):
			if val!=None:
				data_exists=True

		if data_exists==False:
			return None

		# concats the above list of data for any None returns (e.g. [g,r,None] --> [g,r])
		data_arr=[x for x in data_arr if x is not None]	

		# combines all available data into single lists, with 'band' holding the band that this data is from
		ra_arr=[]
		dec_arr=[]
		band=[]

		for element in data_arr:
			ra_arr+=element['data']['ra']
			dec_arr+=element['data']['dec']
			for i in range(0,len(element['data']['ra'])):
				band.append(element['band'])
	
		# set ups an overlay entry, and appends it to the list of entries for this survey
		overlay_data=[]
		for i in range(0,len(ra_arr)):
			overlay_entry={'survey':survey,'position':[ra_arr[i],dec_arr[i]],'radius':None,'corrected':False,'mag':band[i],'marker':marker_type}
			overlay_data.append(overlay_entry)

	return overlay_data

'''
gets data for each survey in 'overlay', and appends the resulting overlay entries into the overal list of overlay detections.
'''
def overlay_selection(image_dict,overlay):
	overlay_arr=[]

	if 'gaia' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='gaia',ra_name='ra',dec_name='dec',mag_name='phot_g_mean_mag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'galex_nuv' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='galex',ra_name='RAJ2000',dec_name='DEJ2000',mag_name='NUVmag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'galex_fuv' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='galex',ra_name='RAJ2000',dec_name='DEJ2000',mag_name='FUVmag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'rosat' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='rosat',ra_name='RAJ2000',dec_name='DEJ2000',mag_name=None,marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'wise' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='wise',ra_name='RAJ2000',dec_name='DEJ2000',mag_name='W1mag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'sdss' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='sdss',ra_name='RA_ICRS',dec_name='DE_ICRS',mag_name='uPmag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'twomass' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='twomass',ra_name='RAJ2000',dec_name='DEJ2000',mag_name='Jmag')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'ztf' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='ztf',ra_name=None,dec_name=None,mag_name=None,overlay_type='tracer',marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'erosita' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='erosita',ra_name='ra',dec_name='dec',mag_name=None,marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'atlas' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='atlas',ra_name=None,dec_name=None,mag_name=None,overlay_type='tracer',marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'gaia_lc' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='gaia_lc',ra_name=None,dec_name=None,mag_name=None,overlay_type='tracer',marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'asassn' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='asassn',ra_name=None,dec_name=None,mag_name=None,overlay_type='tracer',marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data
	if 'crts' in overlay:
		overlay_data=get_overlay(image_dict=image_dict,survey='crts',ra_name=None,dec_name=None,mag_name=None,overlay_type='tracer',marker_type='cross')
		if overlay_data!=None:
			overlay_arr+=overlay_data

	# appends the returned overlay data to the input image_dict
	image_dict['data']['overlay']=overlay_arr

	return image_dict