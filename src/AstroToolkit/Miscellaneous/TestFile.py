from WDPlanetsToolkit.Tools import *

import time
from bokeh.plotting import show

start_time=time.time()

figs=[]

'''
 AR Sco (6050296829033196032,245.4470133276900,-22.8862182469700) is used for all imaging functions as it returns images in all surveys, HU Leo (587316166180416640,141.1853304445800,08.0308343220600)
 is used for everything else.
'''
 
print('testing imaging functions (source):\n')
fig1=getpanstarrsimage(source=6050296829033196032)
fig2=getskymapperimage(source=6050296829033196032)
fig3=getdssimage(source=6050296829033196032)
print('success \n\n\n')

figs.append(fig1)
figs.append(fig2)
figs.append(fig3)

print('testing imaging functions (pos):\n')
fig4=getpanstarrsimage(pos=[245.4470133276900,-22.8862182469700])
fig5=getskymapperimage(pos=[245.4470133276900,-22.8862182469700])
fig6=getdssimage(pos=[245.4470133276900,-22.8862182469700])
print('success \n\n\n')

figs.append(fig4)
figs.append(fig5)
figs.append(fig6)

print('panstarrsquery (source):\n',panstarrsquery(source=587316166180416640),'\n\n\n')
print('skymapperquery (source):\n',skymapperquery(source=587316166180416640),'\n\n\n')
print('gaiaquery: (source):\n',gaiaquery(source=587316166180416640),'\n\n\n')
print('galexquery (source):\n',galexquery(source=587316166180416640),'\n\n\n')
print('rosatquery (source):\n',rosatquery(source=587316166180416640),'\n\n\n')
print('ztfquery (source):\n',ztfquery(source=587316166180416640),'\n\n\n')
print('sdssquery (source):\n',sdssquery(source=587316166180416640),'\n\n\n')
print('wisequery (source):\n',wisequery(source=587316166180416640),'\n\n\n')
print('twomassquery (source):\n',twomassquery(source=587316166180416640),'\n\n\n')

print('panstarrsquery (pos):\n',panstarrsquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('skymapperquery (pos):\n',skymapperquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('gaiaquery (pos):\n',gaiaquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('galexquery (pos):\n',galexquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('rosatquery (pos):\n',rosatquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('ztfquery (pos):\n',ztfquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('sdssquery (pos):\n',sdssquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('wisequery (pos):\n',wisequery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('twomassquery (pos):\n',twomassquery(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')

print('getpanstarrsphot (source):\n',getpanstarrsphot(source=587316166180416640),'\n\n\n')
print('getskymapperphot (source):\n',getskymapperphot(source=587316166180416640),'\n\n\n')
print('getgaiaphot (source):\n',getgaiaphot(source=587316166180416640),'\n\n\n')
print('getgalexphot (source):\n',getgalexphot(source=587316166180416640),'\n\n\n')
print('getsdssphot (source):\n',getsdssphot(source=587316166180416640),'\n\n\n')
print('getwisephot (source):\n',getwisephot(source=587316166180416640),'\n\n\n')
print('gettwomassphot (source):\n',gettwomassphot(source=587316166180416640),'\n\n\n')

print('getpanstarrsphot (pos):\n',getpanstarrsphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('getskymapperphot (pos):\n',getskymapperphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('getgaiaphot (pos):\n',getgaiaphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('getgalexphot (pos):\n',getgalexphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('getsdssphot (pos):\n',getsdssphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('getwisephot (pos):\n',getwisephot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')
print('gettwomassphot (pos):\n',gettwomassphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')

print('getbulkphot (source):\n',getbulkphot(source=587316166180416640),'\n\n\n')

print('getbulkphot (pos):\n',getbulkphot(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')

print('testing ztf lightcurve functions (source):\n')
fig7=getztflc(source=587316166180416640)
print('success \n\n\n')

figs.append(fig7)

print('testing lightcurve functions (pos):\n')
fig8=getztflc(pos=[141.1853304445800,08.0308343220600])
print('success \n\n\n')

figs.append(fig8)

print('testing SED functions (source):\n')
fig9=getsed(source=587316166180416640)
print('success \n\n\n')

figs.append(fig9)

print('testing SED functions (pos):\n')
fig10=getsed(pos=[141.1853304445800,08.0308343220600])
print('success \n\n\n')

figs.append(fig10)

print('getgaiacoords (source):\n',getgaiacoords(source=587316166180416640),'\n\n\n')
print('getgaiasource (pos):\n',getgaiasource(pos=[141.1853304445800,08.0308343220600]),'\n\n\n')

print('testing spectrum functions (source):\n')
fig11=getsdssspectrum(source=587316166180416640)
print('success \n\n\n')

figs.append(fig11)

print('testing spectrum functions (pos):\n')
fig12=getsdssspectrum(pos=[141.1853304445800,08.0308343220600])
print('success \n\n\n')

figs.append(fig12)

print('testing HR diagram function (source):\n')
fig13=gethrd(source=587316166180416640)
print('success \n\n\n')

figs.append(fig13)

print('testing datapage function (source):\n')
fig14=getdatapage(source=587316166180416640)
print('success \n\n\n')

figs.append(fig14)

print('testing datapage function (pos):\n')
fig15=getdatapage(pos=[141.1853304445800,08.0308343220600])
print('success \n\n\n')

figs.append(fig15)

delta_time=round((time.time()-start_time),2)

newline='\n'
print(f'{newline}all tests passed in: {delta_time}s')

for figure in figs:
    show(figure)