from ..Tools import dataquery,getdistance

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from io import StringIO
import math

# this code is adapted from https://stilism.obspm.fr/static/scripts/stilism_dist.py
def getreddening(source=None):
	# get gaia data of source
	gaia_data=dataquery(survey='gaia',source=source)['data']

	gal_lon,gal_lat,parallax,parallax_error=gaia_data['l'][0],gaia_data['b'][0],gaia_data['parallax'][0],gaia_data['parallax_error'][0]
	
	distance=getdistance(parallax)
	distance_plus_err=getdistance(parallax+parallax_error)
	distance_sub_err=getdistance(parallax-parallax_error)

	# perform reddening query
	def query(gal_lon,gal_lat,distance):
		url=f'http://stilism.obspm.fr/reddening?frame=galactic&vlong={gal_lon}&ulong=deg&vlat={gal_lat}&ulat=deg&distance={distance}'
	
		s=requests.Session()
		retries=Retry(total=5,backoff_factor=1)
		s.mount('http://',HTTPAdapter(max_retries=retries))

		r=requests.get(url,allow_redirects=True,timeout=30)

		if r.ok:
			file=StringIO(r.content.decode('utf-8'))
			df=pd.read_csv(file)
			
			return df
		else:
			return None
	
	# get reddening at distance+-error, and use these upper and lower bounds to calculate errors on reddening
	red_df=query(gal_lon=gal_lon,gal_lat=gal_lat,distance=distance)
	red_max_df=query(gal_lon=gal_lon,gal_lat=gal_lat,distance=distance_plus_err)
	red_min_df=query(gal_lon=gal_lon,gal_lat=gal_lat,distance=distance_sub_err)

	if red_df==None or red_max_df==None or red_min_df==None:
		return {'type':'reddening','source':source,'data':None}

	reddening_upper=red_max_df['reddening[mag]'].tolist()[0]-red_df['reddening[mag]'].tolist()[0]
	reddening_lower=red_min_df['reddening[mag]'].tolist()[0]-red_df['reddening[mag]'].tolist()[0]

	# grab values at the central distance
	dist=red_df['distance[pc]'].tolist()[0]
	dist_err=red_df['distance_uncertainty[pc]'].tolist()[0]
	red=red_df['reddening[mag]'].tolist()[0]
	red_err_min=red_df['reddening_uncertainty_min[mag]'].tolist()[0]
	red_err_max=red_df['reddening_uncertainty_max[mag]'].tolist()[0]
	
	# calculate errors
	red_err_lower=math.sqrt(pow(red_err_min,2)+pow(reddening_lower,2))
	red_err_upper=math.sqrt(pow(red_err_max,2)+pow(reddening_upper,2))
		
	reddening_dict={'type':'reddening','source':source,'data':{'dist':distance,'red_dist':dist,'red_dist_err':dist_err,'red':red,'red_upper':round(red+red_err_upper,3),'red_lower':round(red-red_err_lower,3)}}

	return reddening_dict