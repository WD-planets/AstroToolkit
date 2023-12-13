import pandas as pd
import astropy.coordinates as coord
from astroquery.vizier import Vizier
import astropy.units as u
import math
import warnings
from astropy.utils.exceptions import AstropyWarning

catalogue_name='J/A+A/588/A103'

Vizier.ROW_LIMIT = -1

warnings.simplefilter('ignore', category=AstropyWarning)

def ROSATQueryCoords(ObjRa,ObjDec,SearchRadius=3):
	data=[]
	
	#Filter inputs
	if not isinstance(SearchRadius,int) and not isinstance(SearchRadius,float):
		raise Exception('search radius must be a float or an integer')
	
	if not isinstance(ObjRa,float):
		try:
			ObjRa=float(ObjRa)
		except:
			raise Exception('object RA must be a float or an integer')
	if not isinstance(ObjDec,float):
		try:
			ObjDec=float(ObjDec)
		except:
			raise Exception('object DEC must be a float or an integer')
	
	v=Vizier(columns=['**'])
	
	#Sends Vizier query, can return multiple objects in given search radius
	data.append(v.query_region(coord.SkyCoord(ra=ObjRa,dec=ObjDec,unit=(u.deg,u.deg),frame='icrs'),width=SearchRadius*u.arcsec,catalog=catalogue_name))
	if len(data[0])==0:
		print('[ROSAT: ROSATQueryCoords] Error: no ROSAT data found at given coordinates')
		return None
	data=data[0][0].to_pandas()
	data=data.reset_index(drop=True)
	return data
	
def ROSATGetPhotometryCoords(ObjRa,ObjDec,SearchRadius=3):
	data=ROSATQueryCoords(ObjRa,ObjDec,SearchRadius=SearchRadius)
	if not isinstance(data,pd.DataFrame) or data.empty:
		return None
	
	photometry=data[['RAJ2000','DEJ2000','VmagVV10','VTmag','BTmag','VmagBSC','VmagHMXB','VmagLMXB','VmagWD','VsphotWD']].copy()

	return photometry