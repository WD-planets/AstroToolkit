from astroquery.vizier import Vizier
import astropy.coordinates as coord
import pandas as pd
from astropy import units as u

catalogue_name='II/246/out'

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
		print('[TWOMASS: get_data] Error: no 2MASS data found at given coordinates')
		return None
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	return data
	
def get_photometry(ra,dec,radius=3):
	TWOMASSdata=get_data(ra,dec,radius=radius)
	if not isinstance(TWOMASSdata,pd.DataFrame) or TWOMASSdata.empty:
		return None
	
	photometry=TWOMASSdata[['RAJ2000','DEJ2000','_2MASS','Jmag','e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag']].copy()

	return photometry