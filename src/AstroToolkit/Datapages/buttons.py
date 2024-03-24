from bokeh.models import CustomJS, Button
from bokeh.layouts import column

from ..Tools import dataquery
from ..Misc.ProperMotionCorrection import get_adapted_radius

def getinfobuttons(grid_size,source=None,pos=None,simbad_radius=3,vizier_radius=3):
	button_width=round(grid_size/2)
	button_height=round(button_width/3)

	# SIMBAD button
	simbad_button = Button(label="SIMBAD",button_type='primary',height=button_height,width=button_width)	

	if pos!=None:
		ra,dec=pos[0],pos[1]
	elif source!=None:
		gaia_data=dataquery(survey='gaia',source=source)['data']
		ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]
		
		# Scale search radius to include ~26 years of potential proper motion (doesn't actually correct coordinates, just gives a buffer)
		simbad_radius=get_adapted_radius([2016,0],[1990,0],pmra,pmdec,simbad_radius)
		vizier_radius=get_adapted_radius([2016,0],[1990,0],pmra,pmdec,vizier_radius)
	
	# format URL that the button uses
	if pos!=None:
		simbad_url=f'https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={ra}+{dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={simbad_radius}&Radius.unit=arcsec&submit=submit+query&CoordList='
	elif source!=None:
		simbad_url=f'https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={ra}+{dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={simbad_radius}&Radius.unit=arcsec&submit=submit+query&CoordList='
	
	# actual functionality on button click
	simbad_button_js = CustomJS(args=dict(url=simbad_url),code='''
		window.open(url)
	''')
	simbad_button.js_on_event('button_click',simbad_button_js)

	# Vizier button
	vizier_button = Button(label="Vizier",button_type='primary',height=button_height,width=button_width)	
	
	# format URL that the button uses
	if dec>=0:
		vizier_url=f'https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={ra}+{dec}&-c.rs={vizier_radius}&-out.add=_r&-sort=_r&-out.max=$4'
	else:
		vizier_url=f'https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={ra}{dec}&-c.rs={vizier_radius}&-out.add=_r&-sort=_r&-out.max=$4'
	
	# actual functionality on button click
	vizier_button_js = CustomJS(args=dict(url=vizier_url),code='''
		window.open(url)
	''')
	vizier_button.js_on_event('button_click',vizier_button_js)
	
	# Margin = [Top,Right,Bottom,Left]
	simbad_button.margin=[0,round(1/4*grid_size),round(0.025*grid_size),round(1/4*grid_size)]
	vizier_button.margin=[0,round(1/4*grid_size),round(0.025*grid_size),round(1/4*grid_size)]
	
	# combine both buttons into single column object
	buttons=column(simbad_button,vizier_button,align='center')
	buttons.margin=[round(1/4*grid_size),0,round(1/4*grid_size),0]
	
	return buttons