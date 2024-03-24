import matplotlib
import numpy as np
from bokeh.plotting import figure
from bokeh.models import Whisker, ColumnDataSource
from bokeh.transform import linear_cmap
from bokeh.models import CustomJS
from bokeh import events

def plot_lightcurve(data,colour,colours):
	if isinstance(data,dict):
		if data['data']==None:
			print('Note: None object passed to plotlightcurve, suggests no lightcurve data was found in given band.')
			return None
	elif isinstance(data,list):
		available_data=[]
		data_exists=False
		for element in data:
			if element['data']==None:
				print('Note: None object passed to plotlightcurve, suggests no lightcurve data was found in given band.')
			else:
				available_data.append(element)
				data_exists=True
		
		if data_exists==False:
			print('Note: no lightcurve data passed to plotlightcurve, suggests no lightcurve data was found.')
			return None

		data=available_data

	def getcmap(colour):
		if colour=='green':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['greenyellow','forestgreen','greenyellow'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='forestgreen'
		elif colour=='red':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['yellow','red','yellow'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='red'
		elif colour=='blue':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['aqua','royalblue','aqua'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='royalblue'
		elif colour=='black':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['lightgray','black','lightgray'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='black'
		elif colour=='orange':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['gold','orange','gold'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='orange'
		elif colour=='purple':
			colourmap=matplotlib.colors.LinearSegmentedColormap.from_list('',['orchid','darkviolet','orchid'])
			palette=[matplotlib.colors.rgb2hex(c) for c in colourmap(np.linspace(0,1,255))]
			error_colour='darkviolet'
			
		return palette,error_colour

	# make all light curves black if no colours supplied
	if colours==None and isinstance(data,list):
		colours=[]
		for i in range(0,len(data)):
			colours.append('black')
	elif colours!=None and not isinstance(data,list):
		print('Note: colours array will not affect single lightcurves. Use argument colour instead.')

	if isinstance(data,list):
		# strip None elements
		temp=data
		data=[]
		for data_dict in temp:
			if data_dict!=None:
				data.append(data_dict)

		survey_band_list=''
		for data_dict in data:
			survey_band_list+=f"{data_dict['survey']} {data_dict['band']}, "
		survey_band_list=survey_band_list[:-2]

		plot=figure(width=400,height=400,title=f'{survey_band_list} lightcurves',x_axis_label=r'\[\text{Observation Date [days]}\]',y_axis_label=f'{survey_band_list}')
	
	elif isinstance(data,dict):
		survey=data['survey']
		band=data['band']
		plot=figure(width=400,height=400,title=f'{survey} {band} lightcurve',x_axis_label=r'\[\text{Observation Date [days]}\]',y_axis_label=f'{band}')

	def plotting(data_dict,colour,min_time=None):
		if 'hjd_ori' in list(data_dict['data'].keys()):
			time=data_dict['data']['hjd_ori']
		elif 'mjd_ori' in list(data_dict['data'].keys()):
			time=data_dict['data']['mjd_ori']

		# subtract minimum time of combined bands for syncing, otherwise just subtract minimum of current band (i.e. if only plotting one band)
		if min_time==None:
			min_time=min(time)
		for i in range(0,len(time)):
			time[i]=time[i]-min_time

		palette,error_colour=getcmap(colour=colour)			

		mag=data_dict['data']['mag']
		mag_err=data_dict['data']['mag_err']		

		plot.y_range.flipped=True

		# set up data source
		source=ColumnDataSource(data=dict(time=time,mag=mag))

		mapper=linear_cmap(field_name='mag',palette=palette,low=min(mag),high=max(mag))

		# plot points
		plot.circle(x='time',y='mag',source=source,color=mapper,legend_label=f"{data_dict['survey']} {data_dict['band']}")

		# plot error bars
		err_xs=[]
		err_ys=[]
		for x,y,y_err in zip(time,mag,mag_err):
			err_xs.append((x,x))
			err_ys.append((y-y_err,y+y_err))
			
		plot.multi_line(err_xs,err_ys,color=error_colour,legend_label=f"{data_dict['survey']} {data_dict['band']}",level='underlay',line_width=0.5,line_cap='square')

		return plot
	
	# do plotting
	if isinstance(data,list):
		# get total set of times, and subtract minimum for these, so that lightcurves of different bands sync properly
		all_times=[]
		for data_dict in data:
			if 'hjd_ori' in list(data_dict['data'].keys()):
				all_times+=data_dict['data']['hjd_ori']
			elif 'mjd_ori' in list(data_dict['data'].keys()):
				all_times+=data_dict['data']['mjd_ori']

		min_time=min(all_times)

		for data_dict,colour in zip(data,colours):
			plot=plotting(data_dict,colour,min_time)
	else:
		plot=plotting(data,colour)
	
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