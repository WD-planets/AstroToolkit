from bokeh.plotting import figure
from bokeh.models import Label, Range1d
import bokeh

def get_grid(dimensions,plots,grid_size=250):
	# Calculates the area occupied by the input plots (i.e. just normalized by /grid_size)
	unit_area=0
	for i,plot_info_dict in enumerate(plots):
		keys=list(plot_info_dict.keys())

		# checks if the 4 required keys are present in each plot_dict
		if 'name' not in keys or 'width' not in keys or 'height' not in keys or 'figure' not in keys:
			raise Exception("Required keys for plot entries: 'name', 'figure', 'width', 'height'")

		unit_area+=plot_info_dict['width']*plot_info_dict['height']

	dimensions_keys=list(dimensions.keys())
	if 'width' not in dimensions_keys or 'height' not in dimensions_keys:
		raise Exception("Required keys for dimensions: 'width', 'height'")

	grid_width,grid_height=dimensions['width'],dimensions['height']

	# Checks if the unit area calculated above matches the target unit area as given by the dimensions (i.e. all space must be filled)
	if unit_area<grid_width*grid_height:
		raise Exception("Given dimensions must be filled with figures. Pass entries with 'figure':None to fill empty space.")
	elif unit_area>grid_width*grid_height:
		raise Exception('Total area of elements is larger than the given dimensions.')

	returned_plots={}

	#plots is a list of entries in form {'name':name to give plot,'figure':ATKplot,'width':plot_width,'height':plot_height}
	for index,plot_info_dict in enumerate(plots):
		# gets the required parameters
		plot_width=plot_info_dict['width']
		plot_height=plot_info_dict['height']
		plot_name=plot_info_dict['name']
		ATK_plotdict=plot_info_dict['figure']
		
		# if there is an ATK plot dict in the 'figure' key, get its information. If None is in the 'figure' key, this is nonespace used to fill the grid.
		if ATK_plotdict!=None:
			ATK_type=ATK_plotdict['type']
			ATK_plot=ATK_plotdict['plot']

			# if an ATK plot is present and has data, just change its size to fit the grid
			if ATK_plot!=None:
				ATK_plot.height,ATK_plot.width=plot_height*grid_size,plot_width*grid_size
				ATK_plot.sizing_mode='fixed'
			
			# if the actual ATK plot object itself is None (i.e. no data was returned), use this to create a 'missing data' figure)
			else:
				ATK_plot=figure(height=plot_height*grid_size,width=plot_width*grid_size,x_axis_label='placeholder x',y_axis_label='placeholder y')
				ATK_plot.sizing_mode='fixed'
				ATK_plot.x_range,ATK_plot.y_range=Range1d(0,10),Range1d(0,10)

				ATK_plot.xgrid.grid_line_color,ATK_plot.ygrid.grid_line_color,ATK_plot.outline_line_color,ATK_plot.toolbar.logo,ATK_plot.toolbar_location=None,None,None,None,None
				ATK_plot.xaxis.major_label_text_color,ATK_plot.yaxis.major_label_text_color='white','white'
				ATK_plot.xaxis.axis_label_text_color,ATK_plot.yaxis.axis_label_text_color='white','white'

				if ATK_type=='image':
					type_label='Image'
				elif ATK_type=='lightcurve':
					type_label='Light Curve'
				elif ATK_type=='sed':
					type_label='SED'
				elif ATK_type=='spectra':
					type_label='Spectrum'
				elif ATK_type=='powspec':
					type_label='Power Spectrum'
				elif ATK_type=='hrd':
					type_label='HRD'
				
				missing_plot_renderer=Label(x=5,y=5,text=f'Missing {type_label} data',text_align='center',text_font_size='30px')

				ATK_plot.add_layout(missing_plot_renderer)
		else:
			# make a figure and set it to invisible as this is just empty space
			ATK_plot=figure(frame_width=plot_width*grid_size,frame_height=plot_height*grid_size)
			ATK_plot.outline_line_color,ATK_plot.toolbar.logo,ATK_plot.toolbar_location=None,None,None
			ATK_plot.sizing_mode='fixed'

		returned_plots[plot_name]=ATK_plot

	return returned_plots
	