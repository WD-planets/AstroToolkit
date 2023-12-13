import pandas as pd
import astropy.coordinates as coord
from astroquery.vizier import Vizier
import astropy.units as u

catalogue_name='II/335/galex_ais'

Vizier.ROW_LIMIT = -1

def GALEXQueryCoords(ra,dec,radius=3):
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
		print('[GALEX: GALEXQueryCoords] Error: no GALEX data found at given coordinates')
		return None
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	return data
	
def GALEXGetPhotometryCoords(ra,dec,radius=3):
	GALEXdata=GALEXQueryCoords(ra,dec,radius=radius)
	if not isinstance(GALEXdata,pd.DataFrame) or GALEXdata.empty:
		return None
	
	photometry=GALEXdata[['RAJ2000','DEJ2000','objid','NUVmag','e_NUVmag','FUVmag','e_FUVmag']].copy()

	return photometry