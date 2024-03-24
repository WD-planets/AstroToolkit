'''
Changelog:
- added example of config editing to README
- fixed a bug with ATLAS lightcurves very rarely returning empty arrays which led to a crash
- fixed a bug with lightcurve file naming
- improved readfits functionality, can now set default source and ra/dec column names using config or from commandline, and prints similar column names if it couldn't find the given ones.
- added support for CRTS lightcurves to lightcurvequery
- added support for CRTS lightcurves as tracers in imaging overlays
- added note to README about tracers and which ones are 'true' tracers (ztf and crts), and add an example image of wolf 28 being traced by crts
- added crts to readme (lightcurvequery and overlays)
- improved clarity of input description in coord conversion tools in README
- added a getconfigvalue tool that fetches the current value of a given config parameter
- added a getcolumn tool that can be used to grab any column from a .fits file
- added comments throughout main Tools.py source code
- other bug fixes
- plots are now stored as part of a dictionary which stores the plot, as well as basic information such as the survey,source,pos, and type
	- updated showplot, saveplot and export tools so that there is no difference in functionality with the above change
- Now needs BS4, but this should be included already as it is required by astroquery
- updated data dictionaries (returned by dataquery,photquery,spectrumquery, etc.) to be more in line with plot dictionaries returned by plotting tools when returning None
- made spectra legends double-clickable in line with other ATK plots
- improved grid functionality (for datapage creation) in many ways
'''

'''
To-Do:
HIGH PRIORITY
- update README to reflect new changes (plots stored as dictionaries, showplot/saveplot/export now required to use on ATK plot objects, updated grid functionality)

- fix issues with links in github/pypi not working (the ones on github are just me forgetting to format them properly), I think if possible just list features in pypi and have link to github 
  for actual documentation. Not sure if external people (i.e. those outside WDPlanets) can definitely view this repository though?
  
MEDIUM PRIORITY
- add defaults for newly supported surveys to metadata table
- comment new code
- 'plot has no renderer' bokeh warning
- clean up / optimize anything that needs it to make development easier
- make some sort of error logging system that notes any errors encountered during runtime, maybe via a log file of some sort (?)

LOW PRIORITY
- (WEBSITE STILL DOWN?) update ASASA-SN SkyPatrol to v2 (don't like having to github clone it to install, still in beta)
- sort out lightcurve times, i.e. could in theory combine mjd/hjd lightcurves currently (not properly) --> starts to get into combining lightcurves across time domain (Boris spoke about this)
- add setting to config to change default hierarchy in 'any' survey searches (currently images/lightcurves)
'''

# Imports -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

from bokeh.plotting import output_file
import re

from .Misc.file_naming import name_file
from .Misc.input_validation import validateinput

from importlib_resources import files
import configparser
import os

newline='\n'

# Configuration -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

config_file=files('AstroToolkit.Settings').joinpath('config.ini')	

# Create default config file if the file doesn't already exist (i.e. when the package is first installed)
# This handles all default parameters unless they are passed to a tool by the user, in which case these are overwritten
if not os.path.isfile(config_file):
	config=configparser.ConfigParser()

	config.add_section('settings')
	config.set('settings','enable_notifications','True')
	config.set('settings','dataquery_radius','3')
	config.set('settings','photquery_radius','3')
	config.set('settings','bulkphotquery_radius','3')
	config.set('settings','imagequery_size','30')
	config.set('settings','imagequery_overlays','gaia')
	config.set('settings','lightcurvequery_radius','3')
	config.set('settings','atlas_username','None')
	config.set('settings','atlas_password','None')
	config.set('settings','sed_radius','3')
	config.set('settings','spectrum_radius','3')
	config.set('settings','grid_size','250')
	config.set('settings','button_simbad_radius','3')
	config.set('settings','button_vizier_radius','3')
	config.set('settings','plot_size','400')
	config.set('settings','readfits_sourcename','source_id')
	config.set('settings','readfits_coordnames','ra,dec')

	with open(config_file,'w') as configfile:
		config.write(configfile)

# allows user to edit the config file
def editconfig(options):
	edit=configparser.ConfigParser()
	edit.read(config_file)
	settings=edit['settings']		

	accepted_keys=list(settings.keys())

	for key in options:
		if key not in accepted_keys:
			raise Exception(f'Invalid configuration parameter. Accepted parameters are {accepted_keys}.')
	
		settings[key]=str(options[key])
		
	with open(config_file,'w') as configfile:
		edit.write(configfile)
	
# Read config file
def readconfig():
	Config=configparser.ConfigParser()
	Config.read(config_file)

	settings=Config['settings']
	for key in settings:
		if settings[key]=='True':
			settings[key]='1'
		elif settings[key]=='False':
			settings[key]='0'

	config={}		

	# set up dictionary with key:value pairs of the config file
	config['ENABLE_NOTIFICATIONS']=int(settings['enable_notifications'])
	config['DATAQUERY_RADIUS']=float(settings['dataquery_radius'])
	config['PHOTQUERY_RADIUS']=float(settings['photquery_radius'])
	config['BULKPHOTQUERY_RADIUS']=float(settings['bulkphotquery_radius'])
	config['IMAGEQUERY_SIZE']=float(settings['imagequery_size'])
	config['IMAGEQUERY_OVERLAYS']=str(settings['imagequery_overlays'])
	config['LIGHTCURVEQUERY_RADIUS']=float(settings['lightcurvequery_radius'])
	config['ATLAS_USERNAME']=str(settings['atlas_username'])
	config['ATLAS_PASSWORD']=str(settings['atlas_password'])
	config['SED_RADIUS']=float(settings['sed_radius'])
	config['SPECTRUM_RADIUS']=float(settings['spectrum_radius'])
	config['GRID_SIZE']=int(settings['grid_size'])
	config['SIMBAD_RADIUS']=int(settings['button_simbad_radius'])
	config['VIZIER_RADIUS']=int(settings['button_vizier_radius'])
	config['PLOT_SIZE']=int(settings['plot_size'])
	config['READFITS_SOURCENAME']=str(settings['readfits_sourcename'])
	config['READFITS_COORDNAMES']=str(settings['readfits_coordnames']).split(',')
	
	return config

# prints the current value of a given config parameter to the terminal, and returns it
def getconfigvalue(parameter):
	config=readconfig()

	# make parameters returned by readconfig lower case to match the config file
	accepted_parameters=list(config.keys())
	for i,val in enumerate(accepted_parameters):
		accepted_parameters[i]=val.lower()

	if parameter not in accepted_parameters:
		raise Exception(f'{parameter} is not in accepted keys. Accepted keys: {accepted_parameters}')
	else:
		# get current value of parameter
		current_value=config[parameter.upper()]

		print(f'Current Value:{newline}{parameter} : {current_value}')
		return current_value

# Data Query --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def dataquery(survey,pos=None,source=None,radius=None):
	config=readconfig()
	
	# if a radius isn't supplied, use default value from config
	if radius==None:
		radius=config['DATAQUERY_RADIUS']

	# print notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running {survey} dataquery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')

	from .Data.data import survey_map

	validateinput({'survey':survey,'pos':pos,'source':source,'radius':radius},'dataquery')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in dataquery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in dataquery.')

	data=survey_map(survey=survey,pos=pos,source=source,radius=radius)

	return data

# Phot Queries ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def photquery(survey,pos=None,source=None,radius=None):
	config=readconfig()
	
	# get radius from config if one isn't given
	if radius==None:
		radius=config['PHOTQUERY_RADIUS']

	# print notifications of current running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running {survey} photquery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')

	from .Data.photometry import phot_query
	
	validateinput({'survey':survey,'pos':pos,'source':source,'radius':radius},'photquery')
	
	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in photquery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in photquery.')

	photometry=phot_query(survey=survey,pos=pos,source=source,radius=radius)

	return photometry

def bulkphotquery(pos=None,source=None,radius=None):
	config=readconfig()
	
	# get radius from config if one isn't given
	if radius==None:
		radius=config['BULKPHOTQUERY_RADIUS']

	# print notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running bulkphotquery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')

	from .Data.photometry import bulk_query
	
	validateinput({'pos':pos,'source':source,'radius':radius},'bulkphot')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in bulkphotquery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in bulkphotquery.')

	data=bulk_query(pos=pos,source=source,radius=radius)
	
	return data

# Imaging -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def imagequery(survey,pos=None,source=None,size=None,overlays='default',band='g'):
	import numpy as np

	config=readconfig()

	# get default size from config if one isn't supplied
	if size==None:
		size=config['IMAGEQUERY_SIZE']

	# get default overlay list from config if one isn't given
	if overlays=='default':
		overlays=config['IMAGEQUERY_OVERLAYS']
	
	if not isinstance(size,int):
		try:
			size=int(size)
		except:
			size_type=type(size)
			print(f'Invalid size data type. Expected int, got {size_type}.')

	# print notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running {survey} imagequery...{newline}source = {source}{newline}pos = {pos}{newline}size = {size}{newline}overlays = {overlays}{newline}')

	from .Data.imaging import image_correction

	f_return=None
	
	validateinput({'survey':survey,'pos':pos,'source':source},'imagequery')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in imagequery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in imagequery.')

	# set all overlays to enabled
	if overlays=='all':
		overlays='gaia,galex_nuv,galex_fuv,rosat,sdss,twomass,wise,ztf,erosita,atlas,gaia_lc,asassn,crts'

	# split overlay string into a list, and check that they are all valid overlays.
	if overlays!=None:
		overlay_list=overlays.split(',')

		for i in range(0,len(overlay_list)):
			overlay_list[i]=overlay_list[i].lower()
		for i in range(0,len(overlay_list)):
			if overlay_list[i] not in ['gaia','galex_nuv','galex_fuv','rosat','ztf','wise','twomass','sdss','erosita','atlas','gaia_lc','asassn','crts']:
				raise Exception('invalid overlay')
	else:
		overlay_list=[]
	
	# validate some inputs based on the survey. E.g. each survey has a maximum size, and skymapper/panstarrs need to have their 'band' parameters formatted differently
	if survey=='panstarrs':
		if size>1500:
			raise Exception(f'Maximum supported size in {survey} is 1500 arcsec.')
		if not re.match('^[grizy]+$', band):
			raise Exception(f'Invalid {survey} bands. Supported bands are [g,r,i,z,y].')
	
	elif survey=='skymapper':
		if size>600:
			raise Exception(f'Maximum supported size in {survey} is 600 arcsec.')
		if re.match('^[grizuv]+$', band):
			pass
		else:
			raise Exception(f'Invalid {survey} bands. Supported bands are [g,r,i,z,u,v].')
	
		band=list(band)
		temp_string=''
		for i in range(0,len(band)):
			temp_string+=(band[i]+',')
		band=temp_string[:-1]
		
	elif survey=='dss':
		if band!='g':
			print('Note: DSS only supports g band imaging, input band has been ignored.')
		if size>7200:
			print(f'Maximum supported size in {survey} is 7200 arcsec.')

	# do hierarchical query
	if survey=='any':
		image=image_correction(survey='panstarrs',pos=pos,source=source,size=size,band=band,overlay=overlay_list)
		if image['data']==None:
			image=image_correction(survey='skymapper',pos=pos,source=source,size=size,band=band,overlay=overlay_list)
			if image['data']==None:
				image=image_correction(survey='dss',pos=pos,source=source,size=size,band=band,overlay=overlay_list)
				if image['data']==None:
					print('Note: No image found in any supported imaging survey.')
					pass
	
	# do single survey query
	else:
		image=image_correction(survey=survey,pos=pos,source=source,size=size,band=band,overlay=overlay_list)
	
	return image

def plotimage(data):
	config=readconfig()

	# print notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting image...{newline}')	

	from .Plotting.imaging import plot_image
	
	# plot image
	plot=plot_image(image_dict=data)

	# set file name
	if plot!=None:
		filename=name_file(data=data,data_type='ATKimage')
		output_file(filename)

		# scale plot size to values set by config
		plot.width,plot.height=config['PLOT_SIZE'],config['PLOT_SIZE']
	else:
		filename=None

	plot_dict={'type':data['type'],'survey':data['survey'],'source':data['source'],'pos':data['pos'],'plot':plot,'ATKfilename':filename}

	return plot_dict

# HRD Plotting ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def plothrd(source=None,sources=None):
	config=readconfig()

	# prints notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting HRD...{newline}source = {source}{newline}sources = {sources}{newline}')

	from bokeh.plotting import output_file

	from .Plotting.HRD import get_plot
	
	validateinput({'source':source},'plothrd')

	if source!=None and sources!=None:
		raise Exception('Simultaneous source and sources input detected in plothrd.')

	# overlay a single gaia source
	if source!=None:
		plot=get_plot(source=source)
	
	# overlay multiple gaia sources
	elif sources!=None:
		# check data type of sources given (for single source input, this is already done through validateinput above)
		for element in sources:
			if not isinstance(element,int):
				try:
					element=int(element)
				except:
					data_type=type(element)
					raise Exception(f'Incorrect source data type in sources. Expected int, got {data_type}.')
		
		plot=get_plot(sources=sources)
	else:
		raise Exception('source or sources input required in plothrd.')

	if plot!=None:
		# name file
		filename=name_file(data={'source':source,'sources':sources},data_type='ATKhrd')
		output_file(filename)

		# set plot size from config defaults
		plot.width,plot.height=config['PLOT_SIZE'],config['PLOT_SIZE']
	else:
		filename=None

	if source!=None:
		plot_dict={'source':source,'type':'hrd','plot':plot,'ATKfilename':filename}
	elif sources!=None:
		plot_dict={'sources':sources,'type':'hrd','plot':plot,'ATKfilename':filename}

	return plot_dict

# Light Curves ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def lightcurvequery(survey,pos=None,source=None,radius=None,username=None,password=None,sigmaclip=None):
	from .Data.lightcurve_sigma_clip import sigma_clip

	config=readconfig()

	# get radius from config if not given
	if radius==None:
		radius=config['LIGHTCURVEQUERY_RADIUS']
	
	# get username from config if not given (only used in ATLAS queries)
	if username==None:
		username=config['ATLAS_USERNAME']
	# get password from config if not given (only used in ATLAS queries)
	if password==None:
		password=config['ATLAS_PASSWORD']

	# enables notifications of currently running tools if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running {survey} lightcurvequery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')

	from .Data.lightcurves import lightcurve_handling
	
	validateinput({'survey':survey,'pos':pos,'source':source,'radius':radius},'lightcurvequery')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in lightcurvequery')
	elif source==None and pos==None:
		raise Exception('pos or source input required in lightcurvequery.')

	# if any data exists in any band, returns True. Otherwise, returns False
	def check_exists(data):
		data_exists=False
		for element in data:
			if element['data']!=None:
				data_exists=True
		
		return data_exists

	# uses above function to perform heirarchical query if asked for
	if survey=='any':
		data=lightcurvequery(survey='ztf',pos=pos,source=source,radius=radius)
		if check_exists(data)==False:
			data=lightcurvequery(survey='crts',pos=pos,source=source,radius=radius)
			if check_exists(data)==False:
				data=lightcurvequery(survey='atlas',pos=pos,source=source,radius=radius,username=username,password=password)
				if check_exists(data)==False:
					data=lightcurvequery(survey='asassn',pos=pos,source=source,radius=radius)
					if check_exists(data)==False:
						data=lightcurvequery(survey='gaia',pos=pos,source=source,radius=radius)

	# otherwise	perform single survey query
	else:
		data=lightcurve_handling(survey=survey,pos=pos,source=source,radius=radius,username=username,password=password)
	
	# performs sigma clip on data to given sigma level if wanted
	if sigmaclip!=None:
		clipped_data=[]
		for element in data:
			if element['data']!=None:
				clipped=sigma_clip(data=element,sigma=sigmaclip)
				clipped_data.append(clipped)
			else:
				clipped_data.append(None)

		data=clipped_data

	return data

def plotlightcurve(data,colour='black',colours=None):
	config=readconfig()

	# enables notifications of currently running tools if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting lightcurve...{newline}')	

	from .Plotting.lightcurves import plot_lightcurve
	
	if colour not in ['green','red','blue','purple','orange','black']:
		raise Exception('Unsupported colour in plotlightcurve.')

	plot=plot_lightcurve(data=data,colour=colour,colours=colours)

	if plot!=None:
		# names file
		filename=name_file(data=data,data_type='ATKlightcurve')
		output_file(filename)
	
		# sets plot size according to config defaults
		plot.width,plot.height=config['PLOT_SIZE']*2,config['PLOT_SIZE']
	else:
		filename=None
	
	plot_dict={'type':'lightcurve'}

	if isinstance(data,list):
		for i,val in enumerate(data):
			if val!=None:
				plot_dict['survey']=val['survey']
				plot_dict['source']=val['source']
				plot_dict['pos']=val['pos']
	elif isinstance(data,dict):
		plot_dict['survey']=data['survey']
		plot_dict['source']=data['source']
		plot_dict['pos']=data['pos']
	else:
		raise Exception('Expected ATK lightcurve dict or list of ATK lightcurve dicts')

	plot_dict['plot']=plot
	plot_dict['ATKfilename']=filename

	return plot_dict

# SEDs --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def sedquery(pos=None,source=None,radius=None):
	config=readconfig()

	# sets radius from default in config if one isn't supplied
	if radius==None:
		radius=config['SED_RADIUS']	

	# enables notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running sedquery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')	

	from .Data.sed import get_data
	
	data=get_data(pos=pos,source=source,radius=radius)
	
	validateinput({'source':source,'pos':pos},'sedquery')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in dataquery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in dataquery.')

	return data

def plotsed(data):
	config=readconfig()

	# enables notifications of currently running tools if enabled in the config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting SED...{newline}')	

	from .Plotting.sed import plot_sed
	
	plot=plot_sed(sed_data=data)

	if plot!=None:
		# names file
		filename=name_file(data=data,data_type='ATKsed')
		output_file(filename)
	
		# scales plot according to defaults in the config
		plot.width,plot.height=config['PLOT_SIZE']*2,config['PLOT_SIZE']
	else:
		filename=None

	plot_dict={'type':data['type'],'source':data['source'],'pos':data['pos'],'plot':plot,'ATKfilename':filename}

	return plot_dict

# Spectra -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def spectrumquery(survey=None,pos=None,source=None,radius=None):
	config=readconfig()

	# set radius to default in config if one isn't given
	if radius==None:
		radius=config['SPECTRUM_RADIUS']	

	# enables notifications of current running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running {survey} spectrumquery...{newline}source = {source}{newline}pos = {pos}{newline}radius = {radius}{newline}')	

	from .Data.spectra import survey_map
	
	validateinput({'survey':survey,'pos':pos,'source':source,'radius':radius},'spectrumquery')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in spectrumquery.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in spectrumquery.')

	data=survey_map(survey=survey,pos=pos,source=source,radius=radius)
	
	return data

def plotspectrum(data):
	config=readconfig()

	# enables notifications of the currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting spectrum...{newline}')	

	from .Plotting.spectra import get_plot
	
	plot=get_plot(spectrum_dict=data)

	if plot!=None:
		# names file
		filename=name_file(data=data,data_type='ATKspectrum')
		output_file(filename)

		# sets plot sizes to defaults from config
		plot.width,plot.height=config['PLOT_SIZE']*2,config['PLOT_SIZE']
	else:
		filename=None

	plot_dict={'type':'spectra','survey':data['survey'],'source':data['source'],'pos':data['pos'],'plot':plot,'ATKfilename':filename}

	return plot_dict

# Timeseries --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def plotpowspec(data):
	import copy

	config=readconfig()

	# enables notifications of currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Plotting power spectrum...{newline}')	

	from .Plotting.powspec import get_plot
	
	# stops data from being flattened by the powspec plotting
	data_copy=copy.deepcopy(data)

	plot=get_plot(dataset=data_copy)

	if plot!=None:
		# name file
		filename=name_file(data=data,data_type='ATKpowspec')
		output_file(filename)
	
		# set size of file according to defaults in config
		plot.width,plot.height=config['PLOT_SIZE']*2,config['PLOT_SIZE']
	else:
		filename=None

	plot_dict={'type':'powspec'}
	for i,val in enumerate(data):
		if val!=None:
			plot_dict['survey']=val['survey']
			plot_dict['source']=val['source']
			plot_dict['pos']=val['pos']

	plot_dict['plot']=plot
	plot_dict['ATKfilename']=filename

	return plot_dict

def tsanalysis(data):
	config=readconfig()

	# enables notifications of currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Running time series analysis...{newline}')	

	from .Timeseries.ztfanalysis import get_analysis
	
	# load timeseries analysis tool
	get_analysis(dataset=data)

# Data Pages --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def gridsetup(dimensions,plots,grid_size=None):
	config=readconfig()
	
	# sets the grid size to the default given in config if one isn't supplied
	if grid_size==None:
		grid_size=config['GRID_SIZE']

	# enables notifications of currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Creating grid...{newline}width = {dimensions[0]}{newline}height = {dimensions[1]}{newline}grid_size = {grid_size}{newline}')

	from .Datapages.grid import get_grid
	
	plots=get_grid(dimensions=dimensions,plots=plots,grid_size=grid_size)
	
	return plots

def getbuttons(grid_size=None,source=None,pos=None,simbad_radius=None,vizier_radius=None):
	config=readconfig()
	
	# gets default grid size from config if one isn't given
	if grid_size==None:
		grid_size=config['GRID_SIZE']
	# gets default radius to use in SIMBAD query for SIMBAD button if one isn't given
	if simbad_radius==None:
		simbad_radius=config['SIMBAD_RADIUS']
	# gets default radius to use in Vizier query for Vizier button if one isn't given
	if vizier_radius==None:
		vizier_radius=config['VIZIER_RADIUS']

	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Generating info buttons...{newline}source = {source}{newline}pos = {pos}{newline} simbad_radius = {simbad_radius}{newline}vizier_radius = {vizier_radius}{newline}')

	from .Datapages.buttons import getinfobuttons	

	validateinput({'source':source,'pos':pos},'databuttons')

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in getbuttons.')
	elif source==None and pos==None:
		raise Exception('pos or source input detected in getbuttons.')

	# validates data type of some parameters. Since these are only used here, didn't add to validateinput
	if not (isinstance(simbad_radius,int) or isinstance(simbad_radius,float)):
		data_type=type(simbad_radius)
		raise Exception(f'Incorrect simbad_radius data type. Expected float/int, got {data_type}.')
	if not (isinstance(vizier_radius,int) or isinstance(vizier_radius,float)):
		data_type=type(vizier_radius)
		raise Exception(f'Incorrect vizier_radius data type. Expected float/int, got {data_type}.')
	if not (isinstance(grid_size,int) or isinstance(grid_size,float)):
		data_type=type(grid_size)
		raise Exception(f'Incorrect grid_size data type. Expected float/int, got {data_type}.')

	# get buttons
	plot=getinfobuttons(grid_size=grid_size,source=source,pos=pos,simbad_radius=simbad_radius,vizier_radius=vizier_radius)
	
	plot_dict={'type':'buttons','source':source,'pos':pos,'plot':plot}

	return plot_dict

def getmdtable(metadata,pos=None,source=None):
	config=readconfig()

	# enables notifications of currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Generating metadata table...{newline}source = {source}{newline}pos = {pos}{newline}')	

	from .Datapages.metadata import gettable

	if source!=None and pos!=None:
		raise Exception('Simultaneous pos and source input detected in getmdtable.')
	elif source==None and pos==None:
		raise Exception('pos or source input required in getmdtable.')
	
	plot=gettable(metadata_dict=metadata,pos=pos,source=source)
	
	plot_dict={'type':'metadata','source':source,'pos':pos,'plot':plot}

	return plot_dict

# Reddening query ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getreddening(source=None):
	from .Data.reddening import getreddening
	
	if source==None:
		raise Exception('getreddening requires source input.')

	validateinput({'source':source})

	reddening=getreddening(source=source)
	
	return reddening

# File Handling -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def savedata(data):
	config=readconfig()

	# enables notifications for currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Saving data to local storage...{newline}')	

	from .Misc.file_handling import create_file
	
	fname=create_file(data_copy=data)
	
	return fname
	
def readdata(filename):
	config=readconfig()

	# enables notifications for currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Reading data from local storage...{newline}')	

	from .Misc.file_handling import read_file
	
	if not isinstance(filename,str):
		data_type=type(filename)
		raise Exception(f'Incorrect filename datatype. Expected str, got {data_type}.')

	data=read_file(file_name=filename)
	
	return data

# These are taken from bokeh, it is just annoying to have to import separately from bokeh every time you want a plot
def showplot(plot):
	from bokeh.plotting import show

	try:
		if isinstance(plot,dict):
			figure=plot['plot']
		else:
			figure=plot
		show(figure)
	except:
		raise Exception('Unexpected data type passed to showplot, expected ATK plot dict or bokeh plot object')
	
def saveplot(plot):
	from bokeh.plotting import save

	try:
		if isinstance(plot,dict):
			figure=plot['plot']
		else:
			figure=plot
		save(figure)
	except:
		raise Exception('Unexpected data type passed to saveplot, expected ATK plot dict or bokeh plot object')
	
	fname=save(figure)

	print(f'Saved file to {fname}.')

	return fname

# exports plots to PNG. Have to use bokeh save() as a proxy to get the filename, so give an option to keep/delete the .html file.
def export(plot):
	config=readconfig()

	# enables notifications for currently running tool if enabled in config
	if config['ENABLE_NOTIFICATIONS']==1:
		print(f'Exporting plot to PNG...{newline}')

	from bokeh.io import export_png
	from bokeh.plotting import save

	if isinstance(plot,dict):
		try:
			filename=plot['ATKfilename']
			export_png(plot['plot'],filename=f'{filename[:-5]}.png')
		except:
			raise Exception("Unexpected structure of dict passed to export. Required keys: 'plot','ATKfilename'")
	else:
		try:
			# have to save plot as a proxy to get file name, could then rm this but this would potentially delete fils that were not meant to be deleted.
			filename=save(plot)
			export_png(plot,filename=f'{filename[:-5]}.png')
		except:
			raise Exception('Unexpected plot type, expected ATK dict or bokeh plot object.')

	return None

# Miscellaneous -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def correctpm(inputtime,targettime,ra,dec,pmra,pmdec):
	from .Misc.ProperMotionCorrection import PMCorrection
	
	pos=PMCorrection(input=inputtime,target=targettime,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec)

	return pos

def getdistance(parallax):
	# parallax must be in mas
	distance=1/(parallax*10**-3)
	return distance

def convfromdeg(pos):
	from .Misc.coord_conversion import convert_to_hmsdms

	pos=convert_to_hmsdms(pos)
	
	return pos

def convtodeg(pos):
	from .Misc.coord_conversion import convert_to_deg
	
	pos=convert_to_deg(pos)
	
	return pos

def getsources(file_name,col_name=None):
	from .Misc.read_fits import get_source_list
	
	# gets the default source column name from the config if one is not given
	if col_name==None:
		config=readconfig()
		parameter=config['READFITS_SOURCENAME']
	else:
		parameter=col_name

	sources=get_source_list(file_name,parameter)
	return sources

def getpositions(file_name,col_names=None):
	from .Misc.read_fits import get_pos_list

	# here, col_names are ra,dec. Gets these from the config if none are given.
	if col_names==None:
		config=readconfig()
		parameters=config['READFITS_COORDNAMES']
	else:
		parameters=col_names

	pos_list=get_pos_list(file_name,parameters)
	return pos_list

def getcolumn(file_name,col_name):
	from .Misc.read_fits import get_column

	values_list=get_column(file_name,col_name)
	return values_list

# used for installation process, as there is some manual stuff that needs to be done inside package directory
def getpath():
	path = os.path.dirname(__file__)
	print(path)