from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import pandas as pd
from bokeh.plotting import figure
import astropy

catalogue_name='V/154/sdss16'

Vizier.ROW_LIMIT = -1

def get_spectrum_data(ra,dec,radius):
	pos=coords.SkyCoord(ra,dec,unit='deg')
	radius=radius/3600*u.deg

	data=SDSS.get_spectra(pos,radius)
	if data!=None:
		data=data[0][1]
	else:
		return None, None
	
	spectrum_data=data.data
	
	log_wavelength=spectrum_data['loglam']
	wavelength=10**log_wavelength*u.AA
	flux=spectrum_data['flux']*10**-17*u.Unit('erg cm-2 s-1 AA-1')

	return wavelength,flux

def get_plot(ra,dec,radius=3):
	x,y=get_spectrum_data(ra,dec,radius)
	
	if not isinstance(x,astropy.units.quantity.Quantity) or not isinstance(y,astropy.units.quantity.Quantity):
		return None

	plot=figure(width=400,height=400,title="SDSS Spectrum",x_axis_label=r'\[\lambda\text{ }[\text{AA}]\]',y_axis_label=r"\[\text{flux [erg}\text{ cm }^{-2}\text{ s }^{-1}\text{AA}^{-1}]\]")
	plot.line(x,y,color='black',line_width=1)
	
	return plot

def get_data(ra,dec,radius=3):
	data=[]
	
	if not isinstance(radius,int) and not isinstance(radius,float):
		raise Exception('search radius must be a float or an integer')
	
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
	
	v=Vizier(columns=['**'])

	#Sends Vizier query, can return multiple objects in given search radius
	data.append(v.query_region(coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs'),width=radius*u.arcsec,catalog=catalogue_name))
	if len(data[0])==0:
		print('[SDSS: get_data] Error: no SDSS data found at given coordinates')
		return None
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	return data
	
def get_photometry(ra,dec,radius=3):
	SDSSdata=get_data(ra,dec,radius=radius)
	if not isinstance(SDSSdata,pd.DataFrame) or SDSSdata.empty:
		return None
	
	photometry=SDSSdata[['RA_ICRS','DE_ICRS','objID','uPmag','e_uPmag','gPmag','e_gPmag','rPmag','e_rPmag','iPmag','e_iPmag','zPmag','e_zPmag']].copy()

	return photometry