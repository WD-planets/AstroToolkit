import pandas as pd
import math
from astropy.time import Time

from ..Surveys.Gaia import GaiaQueryCoords
from ..Miscellaneous.ProperMotionCorrection import PMCorrection

#Creates Gaia overlay
def getGaiaOverlay(plot,ra,dec,sizeAS,mjd,border,pmra=None,pmdec=None):
	gaiaTime=[2016,0]
	imageTime=Time(mjd,format='mjd').to_datetime()
	imageTime=[imageTime.year,imageTime.month]
	
	if pmra==None and pmdec==None:
		pass
	else:
		sizeAS=PMCorrection(gaiaTime,imageTime,ra,dec,pmra,pmdec,radius=sizeAS)[2]

	detections_made=False
	nearby=GaiaQueryCoords(ra,dec,sizeAS,'dr3')
	#Check if any data was returned by Gaia search
	if isinstance(nearby,pd.DataFrame):
		dataAvailability=True
		nearbyRA=nearby.loc[:,'ra'].tolist()
		nearbyDEC=nearby.loc[:,'dec'].tolist()
		nearbyMagG=nearby.loc[:,'phot_g_mean_mag'].tolist()
		nearbyPMRA=nearby.loc[:,'pmra'].tolist()
		nearbyPMDEC=nearby.loc[:,'pmdec'].tolist()
		
		#Check if any g magnitudes exist
		magCheck=False
		for i in range(0,len(nearbyMagG)):
			if not math.isnan(nearbyMagG[i]):
				magCheck=True
		
		if magCheck==False:
			return plot, detections_made
		
		#correct detection coordinates and list those that were corrected
		correctedIndices=[]
		for i in range(0,len(nearbyRA)):
			if not math.isnan(nearbyPMRA[i]) and not math.isnan(nearbyPMDEC[i]):
				RA_corrected,DEC_corrected=PMCorrection(gaiaTime,imageTime,nearbyRA[i],nearbyDEC[i],nearbyPMRA[i],nearbyPMDEC[i])
				nearbyRA[i]=RA_corrected
				nearbyDEC[i]=DEC_corrected
				correctedIndices.append(i)
		
		#draw detections
		for i in range(0,len(nearby)):
			radiusMultiplier=(nearbyMagG[i]/20.7)
			radius=(border/50+(border/75)**radiusMultiplier)*0.75
			
			if i in correctedIndices and not math.isnan(nearbyMagG[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='red',fill_color=None,legend_label='Gaia g detections (corrected)')
				detections_made=True
			elif i not in correctedIndices and not math.isnan(nearbyMagG[i]):
				plot.circle(x=nearbyRA[i],y=nearbyDEC[i],radius=radius,line_color='red',fill_color=None,line_dash='dotted',legend_label='Gaia g detections (uncorrected)')
				detections_made=True
	else:		
		print('[PanSTARRS: getGaiaOverlay] Note: no Gaia detections found. This will limit accuracy of other detection overlays.')
		return plot, detections_made
	
	return plot, detections_made