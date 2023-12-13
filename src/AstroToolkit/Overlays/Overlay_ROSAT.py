import pandas as pd

from ..Surveys.ROSAT import ROSATQueryCoords

#Creates ROSAT overlay
def getROSATOverlay(plot,ra,dec,sizeAS):
	nearby=ROSATQueryCoords(ra,dec,sizeAS)
	#checks if ROSAT query was successful

	detections_made=False
	if isinstance(nearby,pd.DataFrame):
		nearbyRA=nearby.loc[:,'RAJ2000'].tolist()
		nearbyDEC=nearby.loc[:,'DEJ2000'].tolist()
		
		#draws detections
		for i in range(0,len(nearby)):
			plot.cross(x=nearbyRA[i],y=nearbyDEC[i],color='orange',size=10,legend_label='ROSAT detections (uncorrected)')
			detections_made=True
	
	return plot, detections_made