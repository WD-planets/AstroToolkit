import numpy as np

def validateinput(parameters,tool=None):
	keys=parameters.keys()
	
	# handles survey input validation
	if 'survey' in keys:
		survey=parameters['survey']
		if not isinstance(survey,str):
			data_type=type(survey)
			raise Exception(f'Incorrect survey data type. Expected string, got {data_type}.')
		
		if tool=='imagequery':
			if survey not in ['panstarrs','skymapper','dss','any']:
				raise Exception(f'Unsupported imagequery survey, supported imagequery surveys are [panstarrs,skymapper,dss,any].')
		elif tool=='dataquery':
			if survey not in ['gaia','galex','rosat','panstarrs','skymapper','sdss','twomass','wise','erosita']:
				raise Exception(f'Unsupported dataquery survey. Supported dataquery surveys are [gaia,galex,rosat,panstarrs,skymapper,sdss,twomass,wise,erosita].')
		elif tool=='photquery':
			if survey not in ['gaia','galex','rosat','panstarrs','skymapper','sdss','twomass','wise']:
				raise Exception(f'Unsupported photquery survey. Supported photquery surveys are [gaia,galex,rosat,panstarrs,skymapper,sdss,twomass,wise].')
		elif tool=='lightcurvequery':
			if survey not in ['ztf','atlas','gaia','asassn','crts','any']:
				raise Exception(f'Unsupported lightcurvequery survey. Supported lightcurvequery surveys are [ztf,atlas,gaia,asassn,crts,any].')
		elif tool=='spectrumquery':
			if survey not in ['sdss']:
				raise Exception(f'Unsupported spectrumquery survey. Supported spectrumquery surveys are [sdss].')
	
	# handles pos input validation
	if 'pos' in keys:
		if parameters['pos']!=None:
			pos=parameters['pos']
			if not isinstance(pos,list):
				data_type=type(survey)
				raise Exception(f'Incorrect pos data type. Expected list, got {data_type}.')
			pos_ra,pos_dec=pos[0],pos[1]
		
			if not (isinstance(pos_ra,float) or isinstance(pos_ra,int)) and (isinstance(pos_dec,float) or isinstance(pos_dec,int)):
				pos_ra_data_type=type(pos_ra)
				pos_dec_data_type=type(pos_dec)
				raise Exception(f'Incorrect data type inside pos list. Expected [float/int,float/int], got [{pos_ra_data_type},{pos_dec_data_type}].')
	
	# handles source input validation
	if 'source' in keys:
		if parameters['source']!=None:
			source=parameters['source']
			if not isinstance(source,int) and not isinstance(source,float) and not isinstance(source,str) and not isinstance(source,np.int64):
				data_type=type(source)
				raise Exception(f'Incorrect source data type. Expected int, got {data_type}.')
	
	# handles radius input validation
	if 'radius' in keys:
		radius=parameters['radius']
		if not (isinstance(radius,int) or isinstance(radius,float)):
			data_type=type(radius)
			raise Exception(f'Incorrect radius data type. Expected float/int, got {data_type}.')