from AstroToolkit.Tools import getpositions
from bokeh.plotting import figure,show

def get_plot(file_name):
    pos_list=getpositions(file_name)
  
    plot=figure(width=400,height=400,title='Spatial Distribution',x_axis_label='ra [deg]',y_axis_label='dec [deg]')
    for i in range(0,len(pos_list)):
        plot.circle(x=pos_list[i][0],y=pos_list[i][1],color='black',legend_label='Sources')

    return plot    

if __name__=='__main__':
    fig=getSD(file_name='pre-ztf.fits')
    show(fig)