from bokeh.plotting import figure
from bokeh.models import Whisker, ColumnDataSource, CustomJS, PrintfTickFormatter
from bokeh import events
import math


def plot_sed(sed_data):
    if sed_data==None:
        print('Note: None object passed to plotsed, suggests no SED data was found.')
        return None

    sed_data=sed_data['data']

    # list colours
    colour_arr=['springgreen','royalblue','gold','aquamarine','deepskyblue','orangered','orange','red','black','grey']
    
    plot=figure(width=400,height=400,title='SED',x_axis_label=r'\[\lambda_{\text{eff}}\text{ }[\text{AA}]\]',y_axis_label=r'\[\text{flux [mJy]}\]',x_axis_type='log',y_axis_type="log")
    
    # get list of surveys with data points
    survey_list=[]
    for entry in sed_data:
        survey=entry['survey']
        if not survey in survey_list:
            survey_list.append(survey)

    # set up array of data points from sed_data
    legend=False
    for survey in survey_list:
        x_arr=[]
        y_arr=[]
        err_arr=[]        
        for i in range(0,len(sed_data)):
            if sed_data[i]['survey']==survey:
                x_arr.append(sed_data[i]['wavelength'])
                y_arr.append(sed_data[i]['flux'])
                err_arr.append(sed_data[i]['rel_err'])
        
        # calculate error coords
        err_xs=[]
        err_ys=[]
        for x,y,y_err in zip(x_arr,y_arr,err_arr):
            err_xs.append((x,x))
            err_ys.append((y-y_err,y+y_err))
        
        current_survey_index=survey_list.index(survey)

        # plot data
        for i,v in enumerate(x_arr):
            if not math.isnan(err_arr[i]):
                plot.circle(x=x_arr[i],y=y_arr[i],color=colour_arr[current_survey_index],legend_label=f'{survey}')
            else:
                plot.cross(x=x_arr[i],y=y_arr[i],color=colour_arr[current_survey_index],legend_label=f'{survey} (Upper Limit)',size=7.5)
           
        # plot errors
        plot.multi_line(err_xs,err_ys,color=colour_arr[current_survey_index],legend_label=f'{survey}',level='underlay',line_width=0.5,line_cap='square')
        
        legend=True

    plot.xaxis[0].formatter = PrintfTickFormatter(format="%0e")

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