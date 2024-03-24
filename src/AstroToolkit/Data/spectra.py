from astropy import coordinates as coords
from astroquery.sdss import SDSS
import pandas as pd

def get_spectrum_data(survey,ra,dec,radius=3,source=None,pos=None):
	from astropy import units as u
	
	if survey=='sdss':
		# sets up a skycoord object used in the search (basically just a pos)
		position=coords.SkyCoord(ra,dec,unit='deg')
		radius=radius/3600*u.deg

		# fetches data
		data=SDSS.get_spectra(position,radius,timeout=120)
		if data!=None:
			data=data[0][1]

			spectrum_data=data.data
	
			# formats the data
			log_wavelength=spectrum_data['loglam']
			wavelength=10**log_wavelength
			flux=spectrum_data['flux']*10**-17
		else:
			print(f'Note: {survey} spectrumquery returned no data.')

	if data!=None:
		spectrum_dict={'survey':survey,'source':source,'pos':pos,'data':{'wavelength':wavelength,'flux':flux}}
	else:
		spectrum_dict={'survey':survey,'source':source,'pos':pos,'data':None}

	return spectrum_dict

'''
handles input survey and accounts for proper motion before sending the request
'''
def survey_map(survey,pos=None,source=None,radius=3):
	from ..Data.data import get_survey_times
	from ..Tools import dataquery
	from ..Misc.ProperMotionCorrection import PMCorrection
		
	survey_times=get_survey_times()

	if pos!=None:
		ra,dec=pos[0],pos[1]
	
	elif source!=None:
		gaia_data=dataquery(survey='gaia',source=source)['data']
		# correct for proper motion to survey time
		if gaia_data!=None:
			ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]
			pos_corrected=PMCorrection(survey_times['gaia'],survey_times[survey],ra,dec,pmra,pmdec)
			ra,dec=pos_corrected[0],pos_corrected[1]
		else:
			raise Exception('Gaia source not found.')

	spectrum_dict=get_spectrum_data(survey=survey,ra=ra,dec=dec,radius=radius,source=source,pos=pos)
	
	if spectrum_dict!=None:
		spectrum_dict['type']='spectra'

	return spectrum_dict