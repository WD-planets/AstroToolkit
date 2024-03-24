import pandas as pd
from bokeh.plotting import figure
from bokeh.models import Span, CheckboxGroup, CustomJS, Label, Range1d
from bokeh.layouts import row,column
from bokeh.models import PrintfTickFormatter
from bokeh import events
import math

def get_plot(spectrum_dict):
	if spectrum_dict['data']==None:
		print('Note: None object passed to plotspectrum, suggests no spectrum was found.')
		return None	

	x=spectrum_dict['data']['wavelength']
	y=spectrum_dict['data']['flux']

	plot=figure(width=400,height=400,title="SDSS Spectrum",x_axis_label=r'\[\lambda\text{ }[\text{AA}]\]',y_axis_label=r"\[\text{flux [erg}\text{ cm }^{-2}\text{ s }^{-1}\text{AA}^{-1}]\]",sizing_mode='stretch_both')
	plot.line(x,y,color='black',line_width=1)

	wavelengths=[]
	
	hydrogen={'Hydrogen':[8503,8454,8598,8665,6562,4862,4340,4101.734],'labels':[None,None,None,None,r'\[H\alpha\]',r'\[H\beta\]',r'\[H\gamma\]',r'\[H\delta\]']}
	wavelengths.append(hydrogen)

	helium={'Helium':[4472,4686,4713,4921,5016,5876,6678],'labels':[None,None,None,None,None,None,None]}
	wavelengths.append(helium)
	
	sodium={'Sodium':[8183,8195],'labels':[None,None]}
	wavelengths.append(sodium)
	
	calcium={'Calcium':[3934,3967,8498,8542,8662],'labels':[None,None,None,None,None]}
	wavelengths.append(calcium)
	
	calciumII={'Calcium II':[3608,3854,4109,4383,4737,5165,5636,6191],'labels':[None,None,None,None,None,None,None,None]}
	wavelengths.append(calciumII)

	colours=['red','blue','purple','orange','green','lime','brown']

	spans=[]
	labels=[]
	annotations=[]

	for wavelength_group in wavelengths:
		group_label,annotation_label=list(wavelength_group.keys())[0],list(wavelength_group.keys())[1]
		labels.append(group_label)
		
		span_group=[]		
		annotation_group=[]
	
		index=wavelengths.index(wavelength_group)

		group_size=len(wavelength_group['labels'])
	
		i=0
		for wavelength,annotation in zip(wavelength_group[group_label],wavelength_group[annotation_label]):					
			annotation_height=max(y)+(i/group_size)*0.3*max(y)

			plot.line(x=[0.5*max(x)],y=[0.5*max(x)],legend_label=group_label,line_dash='dashed',line_color=colours[index],line_alpha=0.7,visible=False)

			span_group.append(Span(dimension='height',location=wavelength,line_color=colours[index],line_alpha=0.5,line_dash='dashed'))
			if annotation!=None:
				annotation_group.append(Label(x=wavelength,y=annotation_height,x_offset=2,text=annotation,text_font_size='10pt'))
			else:
				annotation_group.append(Label(x=wavelength,y=0,y_units='screen',text=''))
				
			i+=1

		spans.append(span_group)
		annotations.append(annotation_group)

		plot.renderers.extend(span_group)
		plot.renderers.extend(annotation_group)

	checkbox=CheckboxGroup(
	labels=labels,
	active=list(range(len(wavelengths))),
	)
	
	callback = CustomJS(args=dict(spans=spans,annotations=annotations,checkbox=checkbox),
	code=
	'''
	for(var i=0; i<spans.length; i++){
		for (var j=0; j<spans[i].length; j++){
			spans[i][j].visible = checkbox.active.includes(i);
			annotations[i][j].visible = checkbox.active.includes(i);
		}
	}
	'''
	)

	# when 'active' status of checkbox is changed, trigger callback
	checkbox.js_on_change('active', callback)

	plot.y_range=Range1d(-0.1*max(y),1.4*max(y))

	plot.yaxis[0].formatter = PrintfTickFormatter(format="%0.1e")

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

	layout=row(plot,checkbox)

	return layout