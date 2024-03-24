from bokeh.models import Range1d, CustomJS
import astropy
import numpy as np
from bokeh.plotting import figure, output_file
from bokeh import events

def plot_image(image_dict):
	if image_dict['data']==None:
		print('Note: None object passed to plotimage, suggests no image was found.')
		return None

	size=image_dict['data']['size']

	survey=image_dict['survey']

	# Create axes and plot image
	plot=figure(width=800,height=800,title=f'{survey} Image ({size}")')
	# Fixes scaling issue
	plot.min_border=75

	plot.xgrid.grid_line_color=None
	plot.ygrid.grid_line_color=None
	
	image_data=image_dict['data']['image_header']

	# Get x and y limit (in pixels) from image header
	xlim=image_data['NAXIS1']
	ylim=image_data['NAXIS2']
	
	# not sure if this is a great fix, but DSS sometimes returns images where ylim = xlim + 1, when they should be the same
	if xlim!=ylim:
		xlim=ylim

	# Get points on axes in world coordinates (deg)
	wcs=image_dict['data']['wcs']
	
	# Get array of pixels in x and y
	x_points=np.arange(start=0,stop=xlim+1,step=1)
	y_points=np.arange(start=0,stop=ylim+1,step=1)

	# Apply transformation to coordinates
	coords=wcs.all_pix2world(x_points,y_points,1)
	x_points,y_points=coords[0],coords[1]
	
	# Get image size in degrees
	x_range=max(x_points)-min(x_points)
	y_range=max(y_points)-min(y_points)

	plot.x_range=Range1d(max(x_points),min(x_points))
	plot.y_range=Range1d(min(y_points),max(y_points))

	plot.image(image=[image_dict['data']['image_data']],x=x_points[0],y=y_points[0],dw=x_range,dh=y_range,palette='Greys256',level="image",origin='bottom_right',anchor='bottom_right')

	# apply overlays	
	location=image_dict['data']['location']
	ra,dec=location[0],location[1]
	plot.cross(x=ra,y=dec,color='black',size=10)

	overlay_data=image_dict['data']['overlay']

	# plot overlay detections
	legend=False
	if not overlay_data==None:
		for i in range(0,len(overlay_data)):
			survey=overlay_data[i]['survey']
			mag=overlay_data[i]['mag']
			
			if survey=='gaia':
				colour='red'
				label='Gaia g'
			elif survey=='galex':
				if mag=='NUVmag':
					colour='blue'
					label='GALEX NUV'
				if mag=='FUVmag':
					colour='purple'
					label='GALEX FUV'
			elif survey=='rosat':
				colour='orange'
				label='ROSAT'
			elif survey=='ztf':
				if mag=='g':
					colour='lime'
				elif mag=='r':
					colour='red'
				elif mag=='i':
					colour='aqua'
				label=f'ztf {mag}'
			elif survey=='sdss':
				colour='lime'
				label='SDSS g'
			elif survey=='wise':
				colour='pink'
				label='WISE W1'
			elif survey=='twomass':
				colour='yellow'
				label='2MASS J'
			elif survey=='erosita':
				colour='red'
				label='eROSITA'
			elif survey=='atlas':
				if mag=='c':
					colour='royalblue'
				elif mag=='o':
					colour='blueviolet'
				elif mag=='i':
					colour='violet'
				label=f'ATLAS {mag}'
			elif survey=='gaia_lc':
				if mag=='g':
					colour='lime'
				elif mag=='bp':
					colour='aqua'
				elif mag=='rp':
					colour='red'
				label=f'gaia_lc {mag}'
			elif survey=='asassn':
				if mag=='g':
					colour='green'
				elif mag=='v':
					colour='purple'
				label=f'asassn {mag}'
			elif survey=='crts':
				colour='pink'
				label=f'crts v'

			corrected=overlay_data[i]['corrected']
			if corrected==True:
				line_style='solid'
				label+=' (corrected)'
			else:
				line_style='dotted'
				label+=' (uncorrected)'

			marker_type=overlay_data[i]['marker']			

			# draw detections
			if marker_type=='cross' and corrected==False:
				plot.cross(x=overlay_data[i]['position'][0],y=overlay_data[i]['position'][1],color=colour,size=10,legend_label=label)
				legend=True
			elif marker_type=='cross' and corrected==True:
				plot.square(x=overlay_data[i]['position'][0],y=overlay_data[i]['position'][1],color=colour,size=10,legend_label=label,fill_color=None)
				legend=True
			elif marker_type=='circle':
				plot.circle(x=overlay_data[i]['position'][0],y=overlay_data[i]['position'][1],radius=overlay_data[i]['radius'],line_width=2,line_color=colour,fill_color=None,line_dash=line_style,legend_label=label)
				legend=True				

	if legend==True:
		plot.legend.click_policy="hide"	

		# Double click to hide legend
		toggle_legend_js = CustomJS(args=dict(leg=plot.legend[0]), code='''
				if (leg.visible) {
					leg.visible = false
					}
				else {
					leg.visible = true
				}
		''')
	
		plot.js_on_event(events.DoubleTap, toggle_legend_js)  

	return plot