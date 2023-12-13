import pandas as pd
import math
from astropy.time import Time

from ..Surveys.Gaia import GaiaQueryCoords
from ..Surveys.GALEX import GALEXQueryCoords
from ..Miscellaneous.ProperMotionCorrection import PMCorrection

#Creates Galex overlay
def getGalexOverlayNUV(plot,ra,dec,sizeAS,mjd,border,pmra,pmdec):
	#checks if GALEX search was successful
	gaiaTime=[2016,0]
	imageTime=Time(mjd,format='mjd').to_datetime()
	imageTime=[imageTime.year,imageTime.month]
	galexTime=[2007,0]
	
	if pmra!=None and pmdec!=None:
		sizeAS=PMCorrection(gaiaTime,imageTime,ra,dec,pmra,pmdec,radius=sizeAS)[2]
	
	nearby=GALEXQueryCoords(ra,dec,sizeAS)
	
	detections_made=False
	if isinstance(nearby,pd.DataFrame):
		nearby_gaia=GaiaQueryCoords(ra,dec,sizeAS)
		
		#Plotting of Galex detections + brightness/radius scaling
		nearbyRA=nearby.loc[:,'RAJ2000'].tolist()
		nearbyDEC=nearby.loc[:,'DEJ2000'].tolist()
		nearbyMagUV=nearby.loc[:,'NUVmag'].tolist()
		
		#checks if any NUV magnitudes exist
		magCheck=False
		for i in range(0,len(nearbyMagUV)):
			if not math.isnan(nearbyMagUV[i]):
				magCheck=True
		
		if magCheck==False:
			return plot, detections_made
			
		#fetches gaia data
		if isinstance(nearby_gaia,pd.DataFrame):
			nearbyRA_gaia=nearby_gaia.loc[:,'ra'].tolist()
			nearbyDEC_gaia=nearby_gaia.loc[:,'dec'].tolist()
			nearbyPMRA=nearby_gaia.loc[:,'pmra'].tolist()
			nearbyPMDEC=nearby_gaia.loc[:,'pmdec'].tolist()
			
			#transforms gaia detections back to GALEX time
			for i in range(0,len(nearbyRA_gaia)):
				if not math.isnan(nearbyPMRA[i]) and not math.isnan(nearbyPMDEC[i]):
					RA_corrected,DEC_corrected=PMCorrection(gaiaTime,galexTime,nearbyRA_gaia[i],nearbyDEC_gaia[i],nearbyPMRA[i],nearbyPMDEC[i])
					nearbyRA_gaia[i]=RA_corrected
					nearbyDEC_gaia[i]=DEC_corrected
			
			#checks if there are any GALEX objects within a radius, if there are, applies gaia pm corrections to matching GALEX detection
			correctedIndices=[]
			for i in range(0,len(nearbyRA)):
				for j in range(0,len(nearbyRA_gaia)):
					delta=math.sqrt((nearbyRA[i]-nearbyRA_gaia[j])**2+(nearbyDEC[i]-nearbyDEC_gaia[j])**2)*3600
					if delta<5:
						if not math.isnan(nearbyPMRA[j]) and not math.isnan(nearbyPMDEC[j]):
							RA_corrected,DEC_corrected=PMCorrection(galexTime,imageTime,nearbyRA[i],nearbyDEC[i],nearbyPMRA[j],nearbyPMDEC[j])
							nearbyRA[i]=RA_corrected
							nearbyDEC[i]=DEC_corrected
							correctedIndices.append(i)
							
		#draw detections
		for i in range(0,len(nearby)):
			radiusMultiplier=(nearbyMagUV[i]/20.7)
			radius=(border/50+(border/75)**radiusMultiplier)*0.75
			
			if i in correctedIndices and not math.isnan(nearbyMagUV[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='blue',fill_color=None,legend_label='GALEX NUV detections (corrected)')
				detections_made=True
			elif i not in correctedIndices and not math.isnan(nearbyMagUV[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='purple',fill_color=None,line_dash='dotted',legend_label='GALEX NUV detections (uncorrected)')
				detections_made=True
		
	return plot, detections_made
	
#Creates Galex overlay
def getGalexOverlayFUV(plot,ra,dec,sizeAS,mjd,border,pmra,pmdec):
	gaiaTime=[2016,0]
	imageTime=Time(mjd,format='mjd').to_datetime()
	imageTime=[imageTime.year,imageTime.month]
	galexTime=[2007,0]

	if pmra!=None and pmdec!=None:
		sizeAS=PMCorrection(gaiaTime,imageTime,ra,dec,pmra,pmdec,radius=sizeAS)[2]

	nearby=GALEXQueryCoords(ra,dec,sizeAS)

	detections_made=False
	if isinstance(nearby,pd.DataFrame):
		nearby_gaia=GaiaQueryCoords(ra,dec,sizeAS)
		
		#Plotting of Galex detections + brightness/radius scaling
		nearbyRA=nearby.loc[:,'RAJ2000'].tolist()
		nearbyDEC=nearby.loc[:,'DEJ2000'].tolist()
		nearbyMagUV=nearby.loc[:,'FUVmag'].tolist()
		
		#checks if any FUV magnitudes exist
		magCheck=False
		for i in range(0,len(nearbyMagUV)):
			if not math.isnan(nearbyMagUV[i]):
				magCheck=True
		
		if magCheck==False:
			return plot, detections_made
		
		#fetches gaia data
		if isinstance(nearby_gaia,pd.DataFrame):
			nearbyRA_gaia=nearby_gaia.loc[:,'ra'].tolist()
			nearbyDEC_gaia=nearby_gaia.loc[:,'dec'].tolist()
			nearbyPMRA=nearby_gaia.loc[:,'pmra'].tolist()
			nearbyPMDEC=nearby_gaia.loc[:,'pmdec'].tolist()
			
			#transforms gaia detections back to GALEX time
			for i in range(0,len(nearbyRA_gaia)):
				if not math.isnan(nearbyPMRA[i]) and not math.isnan(nearbyPMDEC[i]):
					RA_corrected,DEC_corrected=PMCorrection(gaiaTime,galexTime,nearbyRA_gaia[i],nearbyDEC_gaia[i],nearbyPMRA[i],nearbyPMDEC[i])
					nearbyRA_gaia[i]=RA_corrected
					nearbyDEC_gaia[i]=DEC_corrected
			
			#checks if there are any GALEX objects within a radius, if there are, applies gaia pm corrections to matching GALEX detection
			correctedIndices=[]
			for i in range(0,len(nearbyRA)):
				for j in range(0,len(nearbyRA_gaia)):
					delta=math.sqrt((nearbyRA[i]-nearbyRA_gaia[j])**2+(nearbyDEC[i]-nearbyDEC_gaia[j])**2)*3600
					if delta<5:
						if not math.isnan(nearbyPMRA[j]) and not math.isnan(nearbyPMDEC[j]):
							RA_corrected,DEC_corrected=PMCorrection(galexTime,imageTime,nearbyRA[i],nearbyDEC[i],nearbyPMRA[j],nearbyPMDEC[j])
							nearbyRA[i]=RA_corrected
							nearbyDEC[i]=DEC_corrected
							correctedIndices.append(i)
		
		#draw detections
		for i in range(0,len(nearby)):
			radiusMultiplier=(nearbyMagUV[i]/20.7)
			radius=(border/50+(border/75)**radiusMultiplier/2)*0.75
			
			if i in correctedIndices and not math.isnan(nearbyMagUV[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='blue',fill_color=None,legend_label='GALEX FUV detections (corrected)')
				detections_made=True
			elif i not in correctedIndices and not math.isnan(nearbyMagUV[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='purple',fill_color=None,line_dash='dotted',legend_label='GALEX FUV detections (uncorrected)')
				detections_made=True
		
	return plot, detections_made