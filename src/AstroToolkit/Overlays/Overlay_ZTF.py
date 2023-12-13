from ..Surveys.ZTF import getData
import pandas as pd

#Creates ZTF overlay
def getZTFOverlay(plot,ra,dec,sizeAS):
	detections_made=False
	data=getData(ra,dec,sizeAS/2,'gri')
	#Checks if data was returned and if so, splits into bands
	if data!=None:
		gData,rData,iData=data[0],data[1],data[2]
	else:
		return plot, detections_made
	
	#checks which bands have data
	do_g,do_r,do_i=False,False,False
	if isinstance(gData,pd.DataFrame):
		do_g=True
	elif isinstance(rData,pd.DataFrame):
		do_r=True
	elif isinstance(iData,pd.DataFrame):
		do_i=True
	
	#draws detections
	if do_g==True:
		nearbyRA=gData.loc[:,'ra'].tolist()
		nearbyDEC=gData.loc[:,'dec'].tolist()

		for i in range(0,len(gData)):
			plot.cross(x=nearbyRA[i],y=nearbyDEC[i],color='aqua',legend_label='ZTF detections (g) (uncorrected)')
			detections_made=True
	elif do_r==True:
		nearbyRA=rData.loc[:,'ra'].tolist()
		nearbyDEC=rData.loc[:,'dec'].tolist()
		
		for i in range(0,len(gData)):
			plot.cross(x=nearbyRA[i],y=nearbyDEC[i],color='aqua',legend_label='ZTF detections (r) (uncorrected)')
			detections_made=True
	elif do_i==True:
		nearbyRA=iData.loc[:,'ra'].tolist()
		nearbyDEC=iData.loc[:,'dec'].tolist()
		
		for i in range(0,len(gData)):
			plot.cross(x=nearbyRA[i],y=nearbyDEC[i],color='aqua',legend_label='ZTF detections (i) (uncorrected)')
			detections_made=True
	
	return plot, detections_made