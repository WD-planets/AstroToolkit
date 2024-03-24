from ..Data.data import survey_map

'''
returns all supported photometry surveys (those used in bulkphot queries)
'''
def get_bulkphot_surveys():
	surveys=['gaia','galex','rosat','sdss','twomass','wise','panstarrs','skymapper']
	return surveys

'''
get photometry for a given survey. uses dataquery, and just filters the resulting data to only include photometry columns
'''
def phot_query(survey,pos=None,source=None,radius=3):
	data_dict=survey_map(survey=survey,pos=pos,source=source,radius=radius)

	skip=False
	if data_dict['data']==None:
		photometry=None
		skip=True

	data=data_dict['data']

	if skip==False:
		if survey=='gaia':
			photometry={key: data[key] for key in ['ra','dec','source_id','phot_g_mean_mag','phot_g_mean_mag_error','phot_bp_mean_mag','phot_bp_mean_mag_error','phot_rp_mean_mag','phot_rp_mean_mag_error']}
		elif survey=='galex':
			photometry={key: data[key] for key in ['RAJ2000','DEJ2000','objid','NUVmag','e_NUVmag','FUVmag','e_FUVmag']}
		elif survey=='rosat':
			photometry={key: data[key] for key in ['RAJ2000','DEJ2000','VmagVV10','VTmag','BTmag','VmagBSC','VmagHMXB','VmagLMXB','VmagWD','VsphotWD']}
		elif survey=='sdss':
			photometry={key: data[key] for key in ['RA_ICRS','DE_ICRS','objID','uPmag','e_uPmag','gPmag','e_gPmag','rPmag','e_rPmag','iPmag','e_iPmag','zPmag','e_zPmag']}
		elif survey=='twomass':
			photometry={key: data[key] for key in ['RAJ2000','DEJ2000','_2MASS','Jmag','e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag']}
		elif survey=='wise':
			photometry={key: data[key] for key in ['RAJ2000','DEJ2000','WISE','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag']}
		elif survey=='panstarrs':
			photometry={key: data[key] for key in ['ra','dec','objID','gMeanPSFMag','gMeanPSFMagErr','rMeanPSFMag','rMeanPSFMagErr','iMeanPSFMag','iMeanPSFMagErr','zMeanPSFMag','zMeanPSFMagErr','yMeanPSFMag','yMeanPSFMagErr']}
		elif survey=='skymapper':
			photometry={key: data[key] for key in ['raj2000','dej2000','object_id','g_psf','e_g_psf','r_psf','e_r_psf','i_psf','e_i_psf','z_psf','e_z_psf','u_psf','e_u_psf','v_psf','e_v_psf']}
		
	phot_dict={'survey':survey,'type':'phot','source':source,'pos':pos,'data':photometry}

	return phot_dict

'''
gets any photometry available in all supported surveys
'''
def bulk_query(pos=None,source=None,radius=3):
	surveys=get_bulkphot_surveys()
	
	# photometry dict is set up as {survey:data,survey:data,...}
	photometry_dict={}
	for survey in surveys:
		photometry=phot_query(survey=survey,pos=pos,source=source,radius=radius)['data']
		photometry_dict[survey]=photometry
	
	# set up final dict
	bulkphot_dict={'type':'bulkphot','source':source,'pos':pos,'data':photometry_dict}
	
	return bulkphot_dict