import requests
from requests.adapters import HTTPAdapter, Retry
from io import BytesIO
import pandas as pd

from bs4 import BeautifulSoup as bs
from numpy import genfromtxt, extract
from io import StringIO
from astropy.coordinates import SkyCoord
from importlib_resources import files
import os
from datetime import datetime
import time

from ..Data.data import get_survey_times

survey_times=get_survey_times()

'''
Probably need to rewrite all of this section at some point, as it is far longer + more convoluted than it probably needs to be
'''

TIMER=files('AstroToolkit.Settings').joinpath('CRTS_TIMER.txt')	

if not os.path.isfile(TIMER):
	with open(TIMER,'w') as timer:
		current_time=datetime.now().strftime("%H:%M:%S")
		timer.write(current_time)

# Adapted from code given by Keith, as CRTS lightcurve interface is not well-suited to script queries
def get_CRTS_lightcurve(ra,dec,radius):
	with open(TIMER,'r') as timer:
		prev_query_time = timer.readline().split(':')

	# calculate time delta in seconds between now and time of last query
	current_time=datetime.now().strftime('%H:%M:%S').split(':')

	# this should fix the case of e.g. 23.59.59 turning to 00.00.01 the next day, which would then report a huge change in time and would allow two queries to be performed in a short time.
	# since there is no differentiation between days, can have issue with it seeing same time 24 hrs apart, and therefore making query wait 15s, but it doesn't matter enough to warrant fixing
	if current_time[0]==0 and prev_query_time[0]==23:
		delta_seconds=0
	else:
		current_time_seconds=int(current_time[2])+60*int(current_time[1])+3600*int(current_time[0])
		prev_query_time_seconds=int(prev_query_time[2])+60*int(prev_query_time[1])+3600*int(prev_query_time[0])

		delta_seconds=abs(current_time_seconds-prev_query_time_seconds)

	# limit frequency of queries
	if delta_seconds < 15:
		wait_time = 15 - delta_seconds
		print(f'CRTS queries are limited to one per 15s. Waiting {wait_time}s before performing query.')
		time.sleep(wait_time)

	requestcoords=SkyCoord(ra,dec, unit="deg",frame='icrs')
    # Get CRTS data   http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv_new.cgi
    # The old url  = http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi?RA=+str(ra)+"&Dec="+str(dec)+"&Rad="+str(radius)+"&DB=photcat&OUT=csv&SHORT=long&PLOT=no"

	# radius / 60 converts from arcseconds to arcmin
	url=f'http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv_new.cgi?RADec={ra} {dec}&Rad={radius/60}&OUT=csv&SHORT=short&DB=photocat'

	with open(TIMER,'w') as timer:
		current_time=datetime.now().strftime("%H:%M:%S")
		timer.write(current_time)

	try:
		cqresult = requests.get(url)
		cqp = bs(cqresult.content, "html.parser")
		filelink = [link.get('href') for link in cqp.find_all('a')]

		if filelink!=[]:
			filelink = filelink[0]
			y = requests.get(filelink)

			if y.status_code ==200: # good transaction
				get_data=genfromtxt(StringIO(y.content.decode(encoding='utf-8')),delimiter=',',names=True,dtype=None)

				# delete rows which are too far away - i.e. beyond search radius)   
				return pd.DataFrame(extract(requestcoords.separation(SkyCoord(get_data['RA'],get_data['Dec'], unit="deg" ,frame='icrs')).arcsecond<radius,get_data))  
			else:
				print("Note: Experiencing issues with CRTS")
				return [None]

		print('Note: CRTS lightcurvequery returned no data.')
		return [None]

	except:
		raise ConnectionError ("error in CRTS_light_curve ra="+str(ra)+"  dec="+str(dec))

def get_data(survey,ra,dec,radius=3,source=None,pos=None,username=None,password=None,pmra=None,pmdec=None):
	if (username!=None or password!=None) and survey!='atlas':
		print('Note: username and password input only required for ATLAS lightcurve queries. Please update these defaults using editconfig.')

	# convert radius into degrees
	radius=radius/3600
	
	if survey=='ztf':
		# failed return should return list of None's equal to number of possible bands
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'r','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'i','data':None}]		

		service='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves'
		url=f'{service}?POS=CIRCLE {ra} {dec} {radius}&BANDNAME=g,r,i&FORMAT=CSV'
	
		# performs query, have to set up some extra querying parameters as ZTF often has issues.
		s=requests.Session()
		retries=Retry(total=5,backoff_factor=1,status_forcelist=[500,502,503,504])
		s.mount('http://',HTTPAdapter(max_retries=retries))
	
		# high timeout value is used here as ZTF data can take a long time to apply to large images in an overlay (as so much data is returned)
		try:
			r=s.get(url,timeout=180)
		except:
			print(f'Note: Experiencing issues with {survey}.')
			return f_return
	
		if r.status_code!=200:
			print(f'Note: Experiencing issues with {survey}.')
			return f_return

		# convert to pandas array (containing all bands)
		data=pd.read_csv(BytesIO(r.content))
		if len(data)==0:
			print(f'Note: {survey} lightcurvequery returned no data.')
			return f_return
	
		# split into separate bands
		gData=data.loc[data['filtercode']=='zg']
		rData=data.loc[data['filtercode']=='zr']
		iData=data.loc[data['filtercode']=='zi']
	
		gData,rData,iData=gData.reset_index(drop=True),rData.reset_index(drop=True),iData.reset_index(drop=True)
	
		data_arr=[{'g':gData},{'r':rData},{'i':iData}]
	
		# set empty data sets to None in dict (e.g. [dict,dict,None]) if no i data is available. note that the above data_arr is a list of dicts, not a single dictionary
		for i in range(0,len(data_arr)):
			for key in data_arr[i]:
				if data_arr[i][key].empty:
					data_arr[i]=None
	
		# sets up list of lightcurve data dictionaries, with None for missing bands, so resulting list always has len() = 3
		bands_arr=['g','r','i']
		data_dict_arr=[]
		for index,band in enumerate(data_arr):
			if band!=None:
				for key in band:
					data=band[key]
					current_band=key

					mag=data.loc[:,'mag'].tolist()
					hjd=data.loc[:,'hjd'].tolist()
					ra=data.loc[:,'ra'].tolist()
					dec=data.loc[:,'dec'].tolist()
					hjd_min=min(hjd)
					hjd_ori=hjd		
					hjd=[x-hjd_min for x in hjd]
					mag_err=data.loc[:,'magerr'].tolist()

					data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':current_band,'data':{'ra':ra,'dec':dec,'hjd':hjd,'hjd_ori':hjd_ori,'mag':mag,'mag_err':mag_err}}
					data_dict_arr.append(data_dict)
			else:
				data_dict_arr.append({'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':bands_arr[index],'data':None})

	# most of code taken from ATLAS' example code: https://fallingstar-data.com/forcedphot/apiguide/
	elif survey=='atlas':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'o','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'c','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'i','data':None}]

		if username==None or password==None:
			raise Exception('username and password required for ATLAS lightcurve query. You can register at: https://fallingstar-data.com/forcedphot/')

		if pmra==None or pmdec==None:
			do_correction=False
		else:
			do_correction=True

		import io,sys,time,re
		from astropy.time import Time
		
		base_url="https://fallingstar-data.com/forcedphot"
		
		r=requests.post(url=f'{base_url}/api-token-auth/',data={'username':username,'password':password})

		if r.status_code==200:
			token=r.json()['token']
			headers={'Authorization':f'Token {token}','Accept':'application/json'}
			
		else:
			print(f'ERROR {r.status_code}.')
			print(r.json())
			return f_return
			
		task_url=None
		while not task_url:
			with requests.Session() as s:
				if do_correction==True:
					# mjd = 57174 is June 1st 2015, which is noted on the site as the start of ATLAS observations. Use 60096 for testing (~6 months at time of writing)
					r=s.post(f'{base_url}/queue/',headers=headers,data={'ra':ra,'dec':dec,'mjd_min': 57174.,'propermotion_ra':pmra,'propermotion_dec':dec,'radec_epoch_year':2000})		
				else:
					r=s.post(f'{base_url}/queue/',headers=headers,data={'ra':ra,'dec':dec,'mjd_min': 57174.})			

				if r.status_code==201:
					task_url=r.json()['url']
				
				elif r.status_code==429:
					message=r.json()['detail']
					print(f'{r.status_code} {message}')
					t_sec=re.findall(r'Available in (\d+) seconds',message)
					t_min=re.findall(r'Available in (\d+) minutes',message)

					if t_sec:
						waittime=int(t_sec[0])
					elif t_min:
						waittime=int(t_min[0])
					else:
						waittime=10
					
					print(f'Waiting {waittime} seconds.')
					time.sleep(waittime)
				else:
					print(f'ERROR {r.status_code}')
					print(r.json())
					sys.exit()
					
		result_url=None
		while not result_url:
			with requests.Session() as s:
				r=s.get(task_url,headers=headers)
				
				if r.status_code==200:
					if r.json()['finishtimestamp']:
						result_url=r.json()['result_url']
						break
						
					elif r.json()['starttimestamp']:
						print(f"Task is running (started at {r.json()['starttimestamp']})")
					else:
						print('Waiting for job to start. Checking again in 10 seconds...')
					time.sleep(10)
				
				else:
					print(f'ERROR {r.status_code}')
					print(r.json())
					sys.exit()
					
		with requests.Session() as s:
			textdata=s.get(result_url,headers=headers).text
			#s.delete(task_url,headers=headers).json()

		df=pd.read_csv(io.StringIO(textdata.replace('###','')),delim_whitespace=True)

		if df.empty:
			return f_return
		
		# subtract minimum mjd from timeseries data for plotting from 0
		mjd=[]
		mjd_ori=df['MJD'].tolist()
		for i in range(0,len(mjd_ori)):
			mjd.append(mjd_ori[i]-min(mjd_ori))
		
		# set negative magnitude values to their absolute values (see https://fallingstar-data.com/forcedphot/apiguide/)
		mag=df['m'].tolist()
		mag_err=df['dm'].tolist()
		for i in range(0,len(mag)):
			mag[i]=abs(mag[i])
		
		ra_arr=df['RA'].tolist()
		dec_arr=df['Dec'].tolist()
		obs_arr=df['Obs'].tolist()
		chi=df['chi/N'].tolist()

		flux=df['uJy'].tolist()
		flux_err=df['duJy'].tolist()

		# filter out any data with fluxes that are below signal-to-noise (<~3), or that have a poor chi squared, and split it into bands
		data_dict_arr=[]
		for band in ['o','c','i']:
			mag_filtered=[]
			mag_err_filtered=[]
			mjd_filtered=[]
			mjd_ori_filtered=[]
			ra_filtered=[]
			dec_filtered=[]

			# do some basic cleaning on poor data
			for i in range(0,len(obs_arr)):
				if obs_arr[i][-1]==band:
					if flux[i]>3 and chi[i]<100 and flux_err[i]<4000:
						mag_filtered.append(mag[i])
						mag_err_filtered.append(mag_err[i])
						mjd_filtered.append(mjd[i])
						mjd_ori_filtered.append(mjd_ori[i])
						ra_filtered.append(ra_arr[i])
						dec_filtered.append(dec_arr[i])

			if len(mag_filtered)>0:
				mag_range = max(mag_filtered)-min(mag_filtered)
				
				mag_filtered2=[]
				mag_err_filtered2=[]
				mjd_filtered2=[]
				mjd_ori_filtered2=[]
				ra_filtered2=[]
				dec_filtered2=[]

				for i in range(0,len(mag_err_filtered)):
					if mag_err_filtered[i]<0.5*mag_range:
						mag_filtered2.append(mag_filtered[i])
						mag_err_filtered2.append(mag_err_filtered[i])
						mjd_filtered2.append(mjd_filtered[i])
						mjd_ori_filtered2.append(mjd_ori_filtered[i])
						ra_filtered2.append(ra_filtered[i])
						dec_filtered2.append(dec_filtered[i])

				data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':band,'data':{'ra':ra_filtered2,'dec':dec_filtered2,'mjd_ori':mjd_ori_filtered2,'mjd':mjd_filtered2,'mag':mag_filtered2,'mag_err':mag_err_filtered2}}
				data_dict_arr.append(data_dict)
			else:
				data_dict_arr.append({'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':band,'data':None})
				
	elif survey=='gaia':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'bp','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'rp','data':None}]		

		from ..Data.data import survey_map
		import numpy as np
		from astropy.time import Time
		import math
		
		# convert radius back into arcseconds
		radius=radius*3600

		data=survey_map(survey='gaia_lc',pos=pos,source=source,radius=radius)
		
		if not isinstance(data,pd.DataFrame):
			print(f'Note: {survey} lightcurvequery returned no data.')
			return f_return
		
		bands=['g','bp','rp']
		
		data_dict_arr=[]
		for band in bands:
			if band=='g':
				# time, flux, flux_err, mag
				cols=['_tab8_5','_tab8_6','_tab8_7','_tab8_9','RA_ICRS','DE_ICRS']
			elif band=='bp':
				# time, flux, flux_err, mag
				cols=['_tab8_11','_tab8_12','_tab8_13','_tab8_15','RA_ICRS','DE_ICRS']
			elif band=='rp':
				# time, flux, flux_err, mag
				cols=['_tab8_16','_tab8_17','_tab8_18','_tab8_20','RA_ICRS','DE_ICRS']

			time,flux,flux_err,mag,ra,dec=data[cols[0]].tolist(),data[cols[1]].tolist(),data[cols[2]].tolist(),data[cols[3]].tolist(),data[cols[4]].tolist(),data[cols[5]].tolist()
		
			# magnitude error formula taken from https://astronomy.stackexchange.com/questions/38371/how-can-i-calculate-the-uncertainties-in-magnitude-like-the-cds-does
			mag_err=[]
			for f,f_e in zip(flux,flux_err):
				mag_err.append((2.5/np.log(10))*(f_e/f))

			bad_indices=[]
			for t in time:
				if math.isnan(t):
					bad_indices.append(time.index(t))
	
			time_f=[]
			mag_f=[]
			mag_err_f=[]
			flux_f=[]
			flux_err_f=[]
			ra_f=[]
			dec_f=[]

			for i in range(0,len(data)):
				if i not in bad_indices:
					time_f.append(time[i])
					mag_f.append(mag[i])
					mag_err_f.append(mag_err[i])
					flux_f.append(flux[i])
					flux_err_f.append(flux_err[i])
					ra_f.append(ra[i])
					dec_f.append(dec[i])

			mjd_ori=[]
			for t in time_f:
				astropy_time=Time(t+2455197.5,format='jd')
				mjd_ori.append(astropy_time.value)
		
			mjd=[]
			for t in mjd_ori:
				mjd.append(t-min(mjd_ori))

			if len(mag)>0:
				data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':band,'data':{'ra':ra_f,'dec':dec_f,'mjd_ori':mjd_ori,'mjd':mjd,'mag':mag_f,'mag_err':mag_err_f}}
				data_dict_arr.append(data_dict)		
			else:
				data_dict_arr.append({'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':band,'data':None})

	elif survey=='asassn':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None}]

		import json
		
		# convert radius to arcminutes
		radius=radius*60

		url=f'https://asas-sn.osu.edu/photometry.json?action=index&controller=photometry&dec={dec}&epochs_max=&epochs_min=&ra={ra}&radius={radius}&rms_max=&rms_min=&sort_by=raj2000&utf8=%E2%9C%93&vmag_max=&vmag_min='

		s=requests.Session()
		retries=Retry(total=5,backoff_factor=1,status_forcelist=[500,502,503,504])
		s.mount('http://',HTTPAdapter(max_retries=retries))

		try:
			r=s.get(url,timeout=60)
		except:
			print(f'Note: Experiencing issues with {survey}.')
			return f_return
		
		data=json.loads(r.content)

		if data['count']>1:
			print(f'Note: {survey} lighcurvequery returned data for multiple objects.')
			
		if data['count']==0:
			print(f'Note: {survey} lightcurvequery returned no data.')
			return f_return

		v_ra_arr=[]
		v_dec_arr=[]
		v_mag=[]
		v_mag_err=[]
		v_hjd_ori=[]
		
		g_ra_arr=[]
		g_dec_arr=[]
		g_mag=[]
		g_mag_err=[]
		g_hjd_ori=[]

		for i in range(0,len(data['results'])):
			link=data['results'][i]['link']			
	
			ra=data['results'][i]['raj2000']
			dec=data['results'][i]['dej2000']

			complete=False			
			while complete==False:
				r=requests.get(link,timeout=60)
			
				lc_data=json.loads(r.content)
				
				for j in range(0,len(lc_data['results'])):
					current_data=lc_data['results'][j]
					camera=current_data['camera'][-1:]
					if camera in ['a','b','c','d','e','f','g','h']:
						v_mag.append(current_data['mag'])
						v_mag_err.append(current_data['mag_err'])
						v_hjd_ori.append(current_data['hjd'])
						v_ra_arr.append(ra)
						v_dec_arr.append(dec)
					elif camera in ['i','j','k','l','m','n','o','p','q','r','s','t']:
						g_mag.append(current_data['mag'])
						g_mag_err.append(current_data['mag_err'])
						g_hjd_ori.append(current_data['hjd'])
						g_ra_arr.append(ra)
						g_dec_arr.append(dec)

				if lc_data['next']!=None:
					link=lc_data['next']
				else:
					complete=True

		v_hjd=[]
		g_hjd=[]
		for time in v_hjd_ori:
			v_hjd.append(time-min(v_hjd_ori))
		for time in g_hjd_ori:
			g_hjd.append(time-min(g_hjd_ori))

		if len(v_mag)==0:
			v_data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':None}
		else:
			v_ra_arr_filtered=[]
			v_dec_arr_filtered=[]
			v_hjd_ori_filtered=[]
			v_hjd_filtered=[]
			v_mag_filtered=[]
			v_mag_err_filtered=[]

			mag_range=max(v_mag)-min(v_mag)	
			for i in range(0,len(v_mag)):
				if v_mag_err[i]<0.5*mag_range:
					v_ra_arr_filtered.append(v_ra_arr[i])
					v_dec_arr_filtered.append(v_dec_arr[i])
					v_hjd_ori_filtered.append(v_hjd_ori[i])
					v_hjd_filtered.append(v_hjd[i])
					v_mag_filtered.append(v_mag[i])
					v_mag_err_filtered.append(v_mag_err[i])
			
			v_data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':{'ra':v_ra_arr_filtered,'dec':v_dec_arr_filtered,'hjd_ori':v_hjd_ori_filtered,'hjd':v_hjd_filtered,'mag':v_mag_filtered,'mag_err':v_mag_err_filtered}}
		
		if len(g_mag)==0:
			g_data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None}
		else:
			g_ra_arr_filtered=[]
			g_dec_arr_filtered=[]
			g_hjd_ori_filtered=[]
			g_hjd_filtered=[]
			g_mag_filtered=[]
			g_mag_err_filtered=[]

			mag_range=max(g_mag)-min(g_mag)	
			for i in range(0,len(g_mag)):
				if g_mag_err[i]<0.5*mag_range:
					g_ra_arr_filtered.append(g_ra_arr[i])
					g_dec_arr_filtered.append(g_dec_arr[i])
					g_hjd_ori_filtered.append(g_hjd_ori[i])
					g_hjd_filtered.append(g_hjd[i])
					g_mag_filtered.append(g_mag[i])
					g_mag_err_filtered.append(g_mag_err[i])
			
			g_data_dict={'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':{'ra':g_ra_arr_filtered,'dec':g_dec_arr_filtered,'hjd_ori':g_hjd_ori_filtered,'hjd':g_hjd_ori_filtered,'mag':g_mag_filtered,'mag_err':g_mag_err_filtered}}

		data_dict_arr=[v_data_dict,g_data_dict]

	elif survey=='crts':
		radius=radius * 3600 # convert radius back into arcsec
		data=get_CRTS_lightcurve(ra,dec,radius)

		if not isinstance(data,pd.DataFrame):
			data_dict_arr=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':None}]
		else:
			ra_arr=data.loc[:,'RA'].tolist()
			dec_arr=data.loc[:,'Dec'].tolist()
			mjd=data.loc[:,'MJD'].tolist()
			mjd_ori=[]
			mjd_min=min(mjd)
			for val in mjd:
				mjd_ori.append(val-mjd_min)
			mag_arr=data.loc[:,'Mag'].tolist()
			mag_err_arr=data.loc[:,'Magerr'].tolist()

			data_dict_arr=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':{'ra':ra_arr,'dec':dec_arr,'mjd_ori':mjd_ori,'mjd':mjd,'mag':mag_arr,'mag_err':mag_err_arr}}]

	return data_dict_arr


'''
corrects for proper motion to the survey time (if source input is used)
'''
def lightcurve_handling(survey,pos=None,source=None,radius=3,username=None,password=None):
	from ..Tools import dataquery
	from ..Misc.ProperMotionCorrection import PMCorrection

	# different surveys can have a different number of bands, so need to return an empty list with that survey's band count (i.e. here, ztf is [None,None,None] for g,r,i)
	if survey=='ztf':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'r','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'i','data':None}]
	elif survey=='atlas':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'o','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'c','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'i','data':None}]
	elif survey=='gaia':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'bp','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'rp','data':None}]
	elif survey=='asassn':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':None},
				  {'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'g','data':None}]
	elif survey=='crts':
		f_return=[{'type':'lightcurve','source':source,'pos':pos,'survey':survey,'band':'v','data':None}]

	if pos!=None:
		ra,dec=pos[0],pos[1]
		pmra,pmdec=None,None
				
		data=get_data(survey=survey,ra=ra,dec=dec,radius=radius,pos=pos,source=source,username=username,password=password)

	elif source!=None:
		gaia_data=dataquery(survey='gaia',source=source)['data']
		if gaia_data!=None:
			ra,dec,pmra,pmdec,ra_j2000,dec_j2000=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0],gaia_data['ra2000'][0],gaia_data['dec2000'][0]
		else:
			return f_return
			
		# correct for proper motion, if the survey is atlas, use in-built proper motion correction.
		if survey=='atlas':
			data=get_data(survey=survey,ra=ra_j2000,dec=dec_j2000,radius=radius,pos=pos,source=source,username=username,password=password,pmra=pmra,pmdec=pmdec)
		else:
			pos_corrected=PMCorrection(input=survey_times['gaia'],target=survey_times[survey],ra=ra,dec=dec,pmra=pmra,pmdec=pmdec)
			ra,dec=pos_corrected[0],pos_corrected[1]
	
			data=get_data(survey=survey,ra=ra,dec=dec,radius=radius,pos=pos,source=source)

	# Had to add this additional check, for some reason ATLAS could return empty data arrays (very rarely), i.e. data['data']['mag']=[] etc. Not sure if would affect other surveys but may
	# aswell apply to all just in case.
	for i,band_data in enumerate(data):
		if band_data['data']!=None:
			if len(band_data['data']['mag'])==0:
				data[i]=f_return[i]

	return data