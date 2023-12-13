from .Overlay_Gaia import getGaiaOverlay
from .Overlay_GALEX import getGalexOverlayFUV
from .Overlay_GALEX import getGalexOverlayNUV
from .Overlay_ROSAT import getROSATOverlay
from .Overlay_ZTF import getZTFOverlay

def overlaySelection(plot,ra,dec,overlay,mjd,sizeAS,border,pmra=None,pmdec=None):
	plot.cross(x=ra,y=dec,color='lime',size=10)
	
	overlay=['gaia']+overlay
	
	for i in range(0,len(overlay)):
		overlay[i]=overlay[i].lower()
	for i in range(0,len(overlay)):
		if overlay[i] not in ['gaia','galex_nuv','rosat','ztf','galex_fuv']:
			raise Exception('invalid overlay')

	adaptive_correction=True
	if pmra==None and pmdec==None:
		adapative_correction=False

	if 'gaia' in overlay:
		if adaptive_correction==True:
			plot,detections_made=getGaiaOverlay(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border,pmra=pmra,pmdec=pmdec)
		else:
			plot,detections_made=getGaiaOverlay(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border)
	if 'galex_fuv' in overlay:
		if adaptive_correction==True:
			plot,detections_made=getGalexOverlayFUV(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border,pmra=pmra,pmdec=pmdec)
		else:
			plot,detections_made=getGalexOverlayFUV(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border)
	if 'galex_nuv' in overlay:
		if adaptive_correction==True:
			plot,detections_made=getGalexOverlayNUV(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border,pmra=pmra,pmdec=pmdec)
		else:
			plot,detections_made=getGalexOverlayNUV(plot=plot,ra=ra,dec=dec,sizeAS=sizeAS,mjd=mjd,border=border)
	if 'rosat' in overlay:
		plot,detections_made=getROSATOverlay(plot,ra,dec,sizeAS)
	if 'ztf' in overlay:
		plot,detections_made=getZTFOverlay(plot,ra,dec,sizeAS)
	
	return plot, detections_made