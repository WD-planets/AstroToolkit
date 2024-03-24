import numpy as np

def convertpos(pos):
	from astropy.coordinates import Angle
	from astropy import units as u

	# Check for a negative declination, only used to later make sure correct signs are shown in the resulting identifier
	ra,dec=pos[0],pos[1]
	if dec<0:
		negativeDec=True
	else:
		negativeDec=False

	ra=Angle(ra,u.degree)
	dec=Angle(dec,u.degree)
	
	# Do unit conversion from deg --> hms/dms
	ra=ra.hms
	dec=dec.dms
	
	# Create ra_arr and dec_arr containing [H,M,S] and [D,M,S]
	ra_arr=np.array([0,0,0],dtype=float)
	dec_arr=np.array([0,0,0],dtype=float)		

	ra_arr[0]=ra[0]
	ra_arr[1]=ra[1]
	ra_arr[2]=ra[2]

	dec_arr[0]=dec[0]
	dec_arr[1]=dec[1]
	dec_arr[2]=dec[2]

	# Convert all negative values to positive (if declination is positive, this is fixed later), this stops the identifier looking like: J...-D-M-S
	ra_str_arr=[]
	for element in ra_arr:
		if element<0:
			element=element*-1

		# Will only retain the SS part from final iteration (which is the only bit you need, since this is the final remainder at the end)
		ra_remainder=element-int(element)

		element=int(element)
		element=str(element).zfill(2) # Force leading zeros to make each element 2 digits long
		ra_str_arr.append(element)

	# Do the same for dec
	dec_str_arr=[]
	for element in dec_arr:
		if element<0:
			element=element*-1
		
		dec_remainder=element-int(element)
		
		element=int(element)
		element=str(element).zfill(2)
		dec_str_arr.append(element)

	# Format remainder: force 2 decimal places, round to 2 decimal places and remove '0.'
	ra_str_arr[2]+=str('{:.2f}'.format(round(ra_remainder,2))[1:])
	dec_str_arr[2]+=str('{:.2f}'.format(round(dec_remainder,2))[1:])
	
	# Writes final identifier and takes account of negative decs.
	if negativeDec==True:
		objRef=f'J{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}-{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}'
	elif negativeDec==False:
		objRef=f'J{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}+{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}'
	else:
		print('objRef Error')
		return None

	return objRef

def convertsource(source):
	from ..Misc.ProperMotionCorrection import PMCorrection as CorrectPM
	from ..Tools import dataquery	

	# Get Gaia data
	gaia_data=dataquery(survey='gaia',source=source)['data']
	if gaia_data!=None:
		ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]
		pos=[ra,dec]
	else:
		return 'SourceNotFound'

	# Correct proper motion to epoch of 2000 used for coordinates. Gaia coordinates are already in ICRS frame of reference so don't need to do anything here.
	pos=CorrectPM([2016,0],[2000,0],ra,dec,pmra,pmdec)

	# Use convertpos() with the corrected object coordinates
	ObjRef=convertpos(pos=pos)

	return ObjRef

def name_file(data,data_type):
	# names files based on their type
	if data_type=='ATKimage':
		source=data['source']
		pos=data['pos']
		survey=data['survey']
	elif data_type=='ATKhrd':
		source=data['source']
		sources=data['sources']
		if sources!=None:
			file_name=f'MultipleSources_ATKhrd.html'
			return file_name
		survey=None
	elif data_type=='ATKpowspec' or data_type=='ATKlightcurve':
		if isinstance(data,list):
			# just grab data from first non-None band, since sources and positions will all be the same
			for band in data:
				if band['data']!=None:
					source=band['source']
					pos=band['pos']
					survey=band['survey']
		else:
			source=data['source']
			pos=data['pos']
			survey=data['survey']
	elif data_type=='ATKbulkphot' or data_type=='ATKsed':
		pos=data['pos']
		source=data['source']
		survey=None
	else:
		source=data['source']
		pos=data['pos']
		survey=data['survey']

	if survey!=None:
		if source!=None:
			obj_ref=convertsource(source=source)
			file_name=f'{obj_ref}_{source}_{survey}_{data_type}.html'
		elif pos!=None:
			file_name=f'{pos[0]}_{pos[1]}_{survey}_{data_type}.html'
	else:
		if source!=None:
			obj_ref=convertsource(source=source)
			file_name=f'{obj_ref}_{source}_{data_type}.html'
		elif pos!=None:
			file_name=f'{pos[0]}_{pos[1]}_{data_type}.html'
		
	return file_name