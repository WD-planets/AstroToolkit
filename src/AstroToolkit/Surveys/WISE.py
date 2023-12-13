from astroquery.vizier import Vizier
import astropy.coordinates as coord
import pandas as pd
from astropy import units as u

catalogue_name='II/311/wise'

Vizier.ROW_LIMIT = -1

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
		print('[WISE: get_data] Error: no WISE data found at given coordinates')
		return None
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	return data
	
def get_photometry(ra,dec,radius=3):
	WISEdata=get_data(ra,dec,radius=radius)
	if not isinstance(WISEdata,pd.DataFrame) or WISEdata.empty:
		return None
	
	photometry=WISEdata[['RAJ2000','DEJ2000','WISE','W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag']].copy()

	return photometry