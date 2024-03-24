from re import S
import astropy.units as u
import numpy as np
import pandas as pd

from ..Tools import bulkphotquery
from .data import get_survey_times

survey_times=get_survey_times()

'''
gets flux using wavelength and zeropoint of a filter, and magnitude through it
'''
def get_flux(mag,zp,wl):
	# sometimes overflows if a bad mag is passed (e.g. -999 for some surveys)
	try:
		flux=zp*10**(-0.4*mag)
	except:
		return None
	c=2.988*10**18
	fnl=1*10**(-23)
	flux=flux/((fnl*c)/wl**2)*1000
	return flux

'''
creates an array of sed points, with each point being a dict {'survey':,'wavelength':,'flux':,'rel_err':}
'''
def format_data(survey,photometry,filter_wavelengths,mag_names,error_names):
	sed_arr=[]
	# gaia has pre-determined zero points, they need to be calculated first for other surveys
	if survey!='gaia':
		for mag_name,error_name,wavelength_ref in zip(mag_names,error_names,filter_wavelengths):
			mag=photometry[mag_name][0]
			mag_err=photometry[error_name][0]
			wavelength=wavelength_ref

			zp=10**((5*np.log10(wavelength)+2.406)/-2.5)
			flux=get_flux(mag=mag,zp=zp,wl=wavelength)

			skip=False
			if flux==None:
				skip=True

			if skip==False:
				rel_err=flux*mag_err/mag
			
				sed_dict={'survey':survey,'wavelength':wavelength,'flux':flux,'rel_err':rel_err}
				sed_arr.append(sed_dict)
	else:
		zero_points=[2.5e-9,4.11e-9,1.24e-9]
		for mag_name,error_name,wavelength_ref,zp in zip(mag_names,error_names,filter_wavelengths,zero_points):
			mag=photometry[mag_name][0]
			mag_err=photometry[error_name][0]
			wavelength=wavelength_ref

			flux=get_flux(mag=mag,zp=zp,wl=wavelength)

			skip=False
			if flux==None:
				skip=True

			if skip==False:
				rel_err=flux*mag_err/mag
			
				sed_dict={'survey':survey,'wavelength':wavelength,'flux':flux,'rel_err':rel_err}
				sed_arr.append(sed_dict)

	return sed_arr

'''
uses bulk phot to gather photometry for all available surveys, and then grabs sed points for each survey that returns data
'''
def get_data(source=None,pos=None,radius=3):
	bulkphot=bulkphotquery(source=source,pos=pos,radius=radius)

	sed_data=[]
	for key in bulkphot['data']:
		if bulkphot['data'][key]!=None:
			phot=bulkphot['data'][key]
			# order of parameters (e.g. mag_names, error_names, zero_points, wavelengths) must match
			if key=='gaia':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[5850.88,5041.61,7690.74],mag_names=['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag'],error_names=['phot_g_mean_mag_error','phot_rp_mean_mag_error','phot_bp_mean_mag_error'])
				sed_data+=data
			if key=='galex':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[2303.37,1548.85],mag_names=['NUVmag','FUVmag'],error_names=['e_NUVmag','e_FUVmag'])
				sed_data+=data
			if key=='sdss':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[3608.04,4671.78,6141.12,7457.89,8922.78],mag_names=['uPmag','gPmag','rPmag','iPmag','zPmag'],error_names=['e_uPmag','e_gPmag','e_rPmag','e_iPmag','e_zPmag'])
				sed_data+=data
			if key=='twomass':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[12350.00,16620.00,21590.00],mag_names=['Jmag','Hmag','Kmag'],error_names=['e_Jmag','e_Hmag','e_Kmag'])
				sed_data+=data
			if key=='wise':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[33526.00,46028.00,115608.00,220883.00],mag_names=['W1mag','W2mag','W3mag','W4mag'],error_names=['e_W1mag','e_W2mag','e_W3mag','e_W4mag'])
				sed_data+=data
			if key=='panstarrs':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[4810.16,6155.47,7503.03,8668.36,9613.60],mag_names=['gMeanPSFMag','rMeanPSFMag','iMeanPSFMag','zMeanPSFMag','yMeanPSFMag'],error_names=['gMeanPSFMagErr','rMeanPSFMagErr','iMeanPSFMagErr','zMeanPSFMagErr','yMeanPSFMagErr'])
				sed_data+=data
			if key=='skymapper':
				data=format_data(survey=key,photometry=phot,filter_wavelengths=[5016.05,6076.85,6076.85,9120.25,3500.22,3878.68],mag_names=['g_psf','r_psf','i_psf','z_psf','u_psf','v_psf'],error_names=['e_g_psf','e_r_psf','e_i_psf','e_z_psf','e_u_psf','e_v_psf'])
				sed_data+=data
	
	sed_dict={'type':'sed','source':source,'pos':pos,'data':sed_data}

	return sed_dict