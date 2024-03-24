def convert_to_deg(pos):
	ra,dec=pos[0],pos[1]
    
	h,m,s=ra[0],ra[1],ra[2]

	converted_ra=15*h+15*m/60+15*s/3600
    
	d,m,s=dec[0],dec[1],dec[2]
	
	if dec[0]>0:
		converted_dec=d+m/60+s/3600
	else:
		converted_dec=d+-1*m/60+-1*s/3600
    
	converted_pos=[converted_ra,converted_dec]

	return converted_pos
    
# this is basically just the file naming code
def convert_to_hmsdms(pos):
	from astropy.coordinates import Angle
	from astropy import units as u
	import numpy as np

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
		converted_pos=[f'{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}',f'-{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}']
	elif negativeDec==False:
		converted_pos=[f'{ra_str_arr[0]}{ra_str_arr[1]}{ra_str_arr[2]}',f'{dec_str_arr[0]}{dec_str_arr[1]}{dec_str_arr[2]}']
	else:
		print('objRef Error')
		return None
	
	return converted_pos