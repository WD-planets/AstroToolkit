from calendar import c
import math

'''
Takes an input year in format [year,month] and corrects coordinates of an object (using its pmra and pmdec) to correct for its proper motion to a target year in format [year,month]
'''
def PMCorrection(input,target,ra,dec,pmra,pmdec):
	# Filter inputs
	if not isinstance(input[0],int):
		raise Exception('input year must be an integer')
	if not isinstance(input[1],int):
		raise Exception('input month must be an integer')
	if not isinstance(target[0],int):
		raise Exception('target year must be an integer')
	if not isinstance(target[1],int):
		raise Exception('target month must be an integer')
	
	if math.isnan(pmra) or math.isnan(pmdec):
		print("Note: Could not retrieve object's pmra/pmdec, and so coordinates were not corrected.")
		return ra,dec

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
		
	if not isinstance(pmdec,float):
		try:
			pmdec=float(pmdec)
		except:
			raise Exception('object PM_RA must be a float or an integer')
		
	if not isinstance(pmdec,float):
		try:
			pmdec=float(pmdec)
		except:
			raise Exception('object PM_DEC must be a float or an integer')
	
	inputYear,inputMonth=input[0],input[1]
	targetYear,targetMonth=target[0],target[1]
	
	# get change in years and change in months
	yearDelta=targetYear-inputYear
	monthDelta=targetMonth-inputMonth

	# correct for proper motion
	ra+=(yearDelta*pmra/3600000+monthDelta*pmra/43200000)*1/math.cos(dec/360*2*math.pi)
	dec+=yearDelta*pmdec/3600000+monthDelta*pmdec/43200000
	
	return ra,dec
	
def get_adapted_radius(input,target,pmra,pmdec,radius):
	inputYear,inputMonth=input[0],input[1]
	targetYear,targetMonth=target[0],target[1]

	if math.isnan(pmra) or math.isnan(pmdec):
		print("Note: Could not retrieve object's pmra/pmdec, and so radius was not adapted.")
		return radius

	# get change in years and change in months
	yearDelta=targetYear-inputYear
	monthDelta=targetMonth-inputMonth
	
	# get time delta in years
	timeDelta=abs(yearDelta+monthDelta/12)	

	# scale radius with proper motion
	radius=radius+math.sqrt((pmra/1000)**2+(pmdec/1000)**2)*timeDelta
	
	return radius