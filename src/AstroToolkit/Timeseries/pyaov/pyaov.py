#!/usr/bin/env python

#
# Python wrapper for period search routines
# (C) Alex Schwarzenberg-Czerny, 2011                alex@camk.edu.pl
# Based on the wrapper scheme contributed by Ewald Zietsman <ewald.zietsman@gmail.com>


import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import csv
'''

Time Series Analysis Package (TSA)
by Alex Schwarzenberg-Czerny (C)2011-2021
Based/Expanded on ESO MIDAS/TSA context by the same author
Modified by Keith Inight 6/7/24 to use wrapping rather than f2py
Modified by Keith Inight 11/7/24 to update the page environment variable so that aovwrapper is found
Modified by Keith Inight 20/7/94 to temporarily update the path environment 
                                    variable so that dynamic link libraries used by
                                    aovwrapper are found
 
Routines For Periodograms/Frequency spectra:
amhw, aovw, atrw, pspw, f_mw, lomw, fgrid
Routines for analysis in Time domain:
covar, fouw, prew 
Auxillary routines:
norm_obs, peak, pldat, plper, simul, totals,

For help type e.g.: help(pyaov.amhw)

NOTE: in all subroutines from times is subtracted EPOCH=(t_min+t_max)/2.

General Reference:
Schwarzenberg-Czerny, A., 1998, Baltic Astronomy, v.7, p.43-69
'''

CTMIN = 5 # global, min. bin population

# a wrapper for the amhw function
def amhw(tin, vin, er, fup, fstep, fr0=0., nh2=3):
    '''
    th,fr,frmax=pyaov.amhw(time, valin, error, fup, fstep, fr0=0., nh2=3)

    Purpose: Returns multiharmonic AOV periodogram, obtained by fitting data
        with a series of trigonometric polynomials. For default nh2=3 this
        is Lomb-Scargle periodogram corrected for constant shift.
    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        nh2[=3]: no. of model parms. (number of harmonics=nh2/2)
        fr0[=0.]: start frequency
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum

    Method:
        General method involving projection onto orthogonal trigonometric 
        polynomials is due to Schwarzenberg-Czerny, 1996. For nh2=2 or 3 it reduces
        Ferraz-Mello (1991), i.e. to Lomb-Scargle periodogram improved by constant 
        shift of values. Advantage of the shift is vividly illustrated by Foster (1995).
    Please quote:
        A. Schwarzenberg-Czerny, 1996, Astrophys. J.,460, L107.   
    Other references:
	G. Foster, 1995, AJ v.109, p.1889 (his Fig.1).
        S. Ferraz-Mello, 1981, AJ v.86, p.619.
	N. R. Lomb, 1976, Ap&SS v.39, p.447.
        A. Schwarzenberg-Czerny, 1998, Baltic Astronomy, 7, 43.
    '''
    nfr = pre("amhw", tin, vin, er, fup, fstep, fr0)
    try:
    
       # th=wrap("amhw",[nh2, 0, fr0, fstep, nfr],[[tin,vin,er]])[3]
        th,_,_,_=aovdrv("AMHW",tin,vin,er,nh2,0,fr0,fstep,nfr,"AOV")
        return post("amhw", th, nfr, fstep, fr0)

    except Exception as e:
        print(e)
        print("amhw: Something unexpected went wrong")
        return 0

# a wrapper for the aovw function
def aovw(tin, vin, er, fup, fstep, fr0=0., nh2=3, ncov=2):
    '''
    th,fr,frmax=pyaov.aovw(time, valin, error, fup, fstep, fr0=0., \
       nh2=3, ncov=2)

    Purpose: Returns AOV periodogram, obtained by phase-folding and binning
        of data.
    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        nh2[=3]:  number of phase bins
        fr0[=0.]: start frequency
        ncov[=2]: number of coverages
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum

    Please quote:
        A. Schwarzenberg-Czerny, 1989, MNRAS 241, 153.
        A. Schwarzenberg-Czerny & J.-Ph. Beaulieu, 2006, MNRAS 365, 165 
    '''

    nfr = pre("aovw", tin, vin, er, fup, fstep, fr0)
    try:
        th,_,_,_=aovdrv("AOVW",tin,vin,er,nh2,ncov,fr0,fstep,nfr,"AOV")
        return post("aovw", th, nfr, fstep, fr0)
    except Exception as e:
        print(e)
        print("aovw: Something unexpected went wrong")
        return 0

# a wrapper for the aovtrw function
def atrw(tin, vin, er, fup, fstep, fr0=0., nh2=30, ncov=2):
    '''
    th,fr,frmax=pyaov.atrw(time, valin, error, fup, fstep, \
         fr0=0., nh2=30, ncov=2)

    Purpose: Returns AOV periodogram for search of transits/eclipses.
        For this purpose data are fitted with a box-like (inverted top-hat) function.

    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        nh2[=30]: number of phase bins
        fr0[=0.]: start frequency
        ncov[=2]: number of bin coverages overlapping in phase
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum

    Please quote:
        A. Schwarzenberg-Czerny, 1989, MNRAS 241, 153.
        A. Schwarzenberg-Czerny & J.-Ph. Beaulieu, 2006, MNRAS 365, 165 
    '''
    
    nfr = pre("atrw", tin, vin, er, fup, fstep, fr0)

    try:

        th,_,_,_=aovdrv("ATRW",tin,vin,er,nh2,ncov,fr0,fstep,nfr,"AOV")
        return post("atrw", th, nfr, fstep, fr0)
    
    
    
    except Exception as e:
        print(e)
        print("atrw: Something unexpected went wrong")
        return 0


def fgrid(time):
    '''
    fup,fstep,fr0=pyaov.fgrid(time)

    Purpose: Evaluate time column and derive a suitable frequency grid
    Input:
        time: numpy array of size (n*1)
    
    Output:  
        fup,fstep,fr0: Frequency maximum, step & minimum, float
    '''

    # check the arrays here, make sure they are all the same size
    try:
        assert time.size >=5
    except AssertionError:
        print('fgrid: Too few time moments')
        return 0

    # maybe something else can go wrong?
    try:
        _,stdout,stderr,parms,return_arrays=wrap("fgrid",[],[time])
        fup,fstep,fr0 = float(parms[0]),float(parms[1]),float(parms[2])
       
        return fup,fstep,fr0

    except Exception as e:
        print(e)
        print("fgrid: Something unexpected went wrong")
        return 0

# a wrapper for the fouw function
def fouw(time, valin, error, frin, nh2=3):
    '''
    frout,dfrout,valout,cof,dcof=pyaov.fouw(time, valin, error, \
                                   frin, bacgnd=0., nh2=3)

    Purpose: fit data with Fourier series, adjusting frequency if needed.
        Returns frequency, pre-whitened residuals,Fourier coefficients and errors
    Input:
        time, valin, error : numpy arrays of size (n*1)
        frin: initial frequency, float
    Optional input:
        nh2[=3]: no. of model parms. (number of harmonics=nh2/2)
    
    Output:  
        frout: final frequency, float
        dfrout: frequency error(if bacgnd correctly set from mhaov|aov|aovtr)
        valout: residuals from fit/data prewhitened with freq
               numpy array of size (n*1)
        cof, dcof: Fourier coefficients & errors,
               numpy arrays of size (m*1) where m = nh2/2*2+2
               cof(m)-float epoch t0, for double use directly time values
               dcof(m)-error of frequency       
    Fourier series:
        fin(ph)=cof(1)+sum_{n=1}^{nh2/2}(cof(2n-1)Cos(ph*n)+cof(2n)*Sin(ph*n))
        where ph=2pi*frout*(t-t0) and t0=(max(t)+min(t))/2

    Please quote:
        A.Schwarzenberg-Czerny, 1995, Astr. & Astroph. Suppl, 110,405
    '''

    # check the arrays here, make sure they are all the same size
    try:
        assert time.size == valin.size == error.size
    except AssertionError:
        print('fouw: Input arrays must be same size')
        return 0

    # maybe something else can go wrong?
    try:

        _,stdout,stderr,parms,return_arrays=wrap("fouw",[frin,nh2],[time, valin, error,])
        frout=parms[0]
        dfrout=parms[1]
        valout=np.array([float(r) for r in return_arrays[0]])
        cof=np.array([float(r) for r in return_arrays[1]])
        dcof=np.array([float(r) for r in return_arrays[2]])
        return frout,dfrout,valout,cof,dcof

    except Exception as e:
        print(e)
        print("fouw: Something unexpected went wrong")
        return 0

def f_mw(tin, vin, er, fup, fstep, fr0=0.):
    '''
    th,fr,frmax=pyaov.f_mw(time, valin, error, fup, fstep, fr0=0.)

    Purpose: Returns generalized L-S periodogram
        of data, using the original Ferraz-Mello (1981) algorithm.
    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        fr0[=0.]: start frequency
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum

    Reference:        
        S. Ferraz-Mello, 1981, Astronomical Journal 86, 619.
        A. Schwarzenberg-Czerny, 1998, Baltic Astronomy, 7, 43.
    '''

    nfr = pre("f_mw", tin, vin, er, fup, fstep, fr0)
    try:
        th,_,_,_=aovdrv("F_MW",tin,vin,er,0,0,fr0,fstep,nfr,"RAW")
        return post("f_mw", th, nfr, fstep, fr0)
    
    except Exception as e:
        print(e)
        print("f_mw: Something unexpected went wrong")
        return 0

def lomw(tin, vin, er, fup, fstep, fr0=0.):
    '''
    th,fr,frmax=pyaov.lomw(time, valin, error, fup, fstep, fr0=0., \
       nh2=3, ncov=2)

    Purpose: Returns AOV periodogram, obtained by phase-folding and binning
        of data.
    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        fr0[=0.]: start frequency
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum

    Reference: 
        N.R. Lomb, 1976, Astrophysics and Space Science, 39, 447. 
        A. Schwarzenberg-Czerny, 1998, Baltic Astronomy, 7, 43.
    '''

    nfr = pre("lomw", tin, vin, er, fup, fstep, fr0)
    try:
        th,_,_,_=aovdrv("LOMW",tin,vin,er,0,0,fr0,fstep,nfr,"RAW")
        return post("lomw", th, nfr, fstep, fr0)
    except Exception as e:
        print(e)
        print("lomw: Something unexpected went wrong")
        return 0






def normalize(x,error,mean=0.,var=1.):
    '''
    y = pyaov.normalize(x,error,mean=0.,var=1.)

    Purpose: converts series to one of specified mean and variance var.
    Input:
        x: input time series, numpy array (n*1)
        error: used to calculate weights w = 1/error**2,
               observations with null error are ignored
    Optional Input:
        mean,var: final mean and variance of the normalized series. 
             For var==0 only mean is shifted.

    Output:  
        x: normalized time series, numpy array (n*1)

    '''
    
    # check the arrays here, make sure they are all the same size
    try:
        assert x.size > 1, 'Normalize:Input arrays length must be >1'
        assert x.size == error.size,'Normalize:Input arrays length must be equal'
        assert var >= 0.,'Normalize:negative target variance'
    except AssertionError:
        return 0

    # maybe something else can go wrong?
    try:
        w = error
        ndx = w <= 0.
        w[ndx] = 1.
        w = 1./np.square(w)
        w[ndx] = 0.
        n = len(w)-len(ndx)
        
        sw = np.sum(w)
        av = np.sum(x*w)/sw
        y = x-av
        vr = np.sum(np.square(y)*w)
        va = vr/sw*n/(n-1)
        if (var>0. and va>0.):
          y=y*np.sqrt(var/va)
        y=y+mean     
        return y, n, sw, av, vr

    except Exception as e:
        print (e)
        print ("Normalize: Something unexpected went wrong!!")
        return 0

def peak(f):
    '''
    xm,fm,dx = pyaov.peak(f)
    from ort
     Scan EVENSAMPLED series for a peak
     Result location & width is in units of f indices (starting from 0)
     INPUT:
       f:- values at even spaced arguments, numpy array (n*1)
     OUTPUT:
       xm-peak location (first element: 0 in Py, 1 in F90)
       fm-peak value
       dx-peak halfwidth at zero level
     SPECIAL CONDITIONS:
       dx<0 no valid peak with positive value
     METHOD:
       Fits parabola to top 3 points
       (C) Alex Schwarzenberg-Czerny, 1999-2005 alex@camk.edu.pl
  '''
    try:
        assert f.size >0
    except AssertionError:
        print('peak: Too few components')
        return 0

    # maybe something else can go wrong?
    try:
        _,stdout,stderr,parms,return_arrays=wrap("peak",[],[f])
        xm, fm, dx = float(parms[0]),float(parms[1]),float(parms[2])

    # NOTE: aov.aovsub.peak follows Fortran indexing convention
    # For Python 1 must be subtracted from xm output
        return xm-1., fm, dx

    except Exception as e:
        print(e)
        print("peak: Something unexpected went wrong")
        return 0
        
    # make a plot of data
def pldat(time, value):
    '''
    pyaov.pldat(time, value)

    Purpose:
        Plots periodogram and data folded with its peak frequency

    Input:
       time, value: numpy arrays of size (n*1)
    
    Output:
        Plot window. Close it to proceed further.
    '''

    pl.figure(figsize=(9,9))
    
    pl.subplot(211)
    pl.plot(time, value, 'k.')
#    pl.title('Data')
    pl.xlabel('Time')
    pl.ylabel('Data Value')

    pl.show()

    # make a periodogram plot
def plper(nam,frmax, time, value, freqs, th):
    '''
    pyaov.plper(frmax, time, value, freqs, th)

    Purpose:
        Plots periodogram and data folded with its peak frequency

    Input:
        frmax: the frequency used for phase folding
           set frmax=0 to avoid phase folded plot
        time, value: numpy arrays of size (n*1)
        freqs, th: frequencies & periodogram (m*1)
    
    Output:
        Plot window. Close it to proceed further.
    '''

    pl.figure(figsize=(9,9))
    
    if frmax>0. :   # plot phase folded data
       pl.subplot(211)
       ph=np.mod(time*frmax,1.)
       pl.plot(ph, value, 'k.')
       pl.xlabel('Phase')
       pl.ylabel('Data Value')
       pl.subplots_adjust( hspace=0.3)

    pl.subplot(212) # plot periodogram
    pl.plot(freqs, th, 'k')
    pl.title(nam)
    pl.ylabel('Statistics')
    pl.xlabel('Frequency')

    pl.show()

def post(nam, th, nfr, fstep, fr0): 
    '''
 Auxillary routine for post-processing of periodograms  
 th, freqs, frmax = post(nam, th, nfr, fstep, fr0)
    '''     
    xm, fm, dx = peak(th)  
    frmax = fr0 + xm * fstep
    print(nam+': Peak of %f at frequency %f' % (fm,frmax))
# make an array that contains the frequencies too
    freqs = np.arange(nfr)*fstep+fr0
    return th, freqs, frmax

def pre(nam, tin, vin, er, fup, fstep, fr0):
# test input for periodograms       
    '''
 Auxillary routine for testing input oto periodogram routines. 
 nfr = pre(nam, tin, vin, er, fup, fstep, fr0)
    '''     
    nobs = tin.size
    assert nobs == vin.size == er.size, \
        nam+': Input arrays must be the same dimensions'
    nfr = int((fup - fr0)/fstep+1.5)
    assert nfr>0, nam+': wrong frequency range'
    if (tin.max()-tin.min())*fstep>0.2:    
        print(nam+': warning: undersampling in frequency')  
    return nfr    

# a wrapper for the prew function
def prew(time, valin, error, frin, nh2=3):
    '''
    frout,dfrout,valout=pyaov.prew(time, valin, error, frin, nh2=3)

    Purpose:
        Improve trig orthogonal polynomial fit by tweaking of frequency.
        Returns frequency and pre-whitened residuals. 
        Note: valin-valout returns a fitted interpolation trig polynomial.

    Input:
        time, valin, error : numpy arrays of size (n*1)
        frin: initial frequency, float
    Optional input:
        nh2[=3]: no. of model parms. (number of harmonics=nh2/2)
    
    Output:  
        frout: final frequency, float
        dfrout: its error
        valout: residuals from fit/data prewhitened with freq
               numpy array of size (n*1)

    Please quote:
        A.Schwarzenberg-Czerny, 1996, Astrophys. J.,460, L107.   
    '''

    # check the arrays here, make sure they are all the same size
    try:
        assert time.size == valin.size == error.size
    except AssertionError:
        print('prew: Input arrays must be the same dimensions')
        return 0

    # maybe something else can go wrong?
    try:
        _,stdout,stderr,parms,return_arrays=wrap("prew",[frin,nh2],[time,valin,error])
    
        th=np.array([float(r) for r in return_arrays[0]])
        return     float(parms[0]),float(parms[1]),th


    except Exception as e:
        print(e)
        print("prew: Something unexpected went wrong")
        return 0
        
# a wrapper for the powspw function
def pspw(tin, vin, er, fup, fstep, fr0=0.):
    '''
    th,fr,frmax=pyaov.pspw(time, valin, error, fup, fstep, fr0=0.)

    Purpose: Returns Power Spectrum. Principal use for calculation of window function.  
        For periodogram we recommend amhv routine with nh2=2 instead.
    Input:
        time, valin, error : numpy arrays of size (n*1)
        fup: frequency to stop calculation at, float
        fstep: size of frequency steps, float
    Optional input:
        fr0[=0.]: start frequency
    
    Output:
        th,fr: periodogram values & frequencies: numpy arrays of size (m*1)
              where m = (fup-fr0)/fstep+1
        frmax: frequency of maximum
    Method:
    	For unit weights exactly corresponds to Deeming DPS. Weights are 
    	assumed to be proportional to number of fictional repeated measurements
    	and so accountend within Deeming scheme.
    Reference:
        T. J. Deeming, 1975, Ap&SS, v36, pp.137-158
        A. Schwarzenberg-Czerny, 1998, Baltic Astronomy, 7, 43.

    '''
    nfr = pre("pspw", tin, vin, er, fup, fstep, fr0)
    try:
        th,_,_,_=aovdrv("PSPW",tin,vin,er,0,0,fr0,fstep,nfr,"RAW")
        return post("pspw", th, nfr, fstep, fr0)
    except Exception as e:
        print(e)
        print("pspw: Something unexpected went wrong")
        return 0
        
def refine(xd, yd, ed, fr): # 
    '''
    fr,dfr,valout,nh2 = refine(xd, yd, ed, fr)
    Refine frequency by orthogonal polynomials fit initial data from box text
    '''
    try:
        n=xd.size
        assert n >0 and n == yd.size and n == ed.size, \
           "refine: wrong array sizes"   
        assert fr > 0., "refine: frequency must be positive"
    except AssertionError:
        return 0
        
    try:
        nh2new = 1
        ndx = ed > 0.
        n = len(ndx)
        snew = np.sum(np.square(yd[ndx]/ed[ndx]))
        next = True
        valout = np.array(yd)
        frout = (fr)

        while next: # fit a longer orthogonal polynomial series 
           nh2old=(nh2new)
           sold=(snew)
           nh2new += 2
            
           fr, dfr, val=prew(xd, yd, ed, frout, nh2=nh2new)
              
           snew = np.sum(np.square(val[ndx]/ed[ndx]))
           fnew = (sold-snew)*(n-nh2new)/(2*snew)
#     stop if (S_{j-1} - S_j)/S_j < F(1,N-j-1) = Fisher(N-j-1)
#     where S=sum((O-C)^2*w)
           if ( (fnew > Fisher(n-nh2new)) or (nh2new < 4) ): # still improving fit
              next = True
              dfrout = (dfr)
              frout = (fr)
              valout=np.array(val)
           else:
              next = False
        return frout,dfrout,valout,nh2new-2 
    except Exception as e:
        print(e)
        print("pspw: Something unexpected went wrong")
        return 0
        
def Fisher(x):
    # returns Fisher-Snedecor test value, F(1,x) for q=0.95           
    # holds for integer x, 2<= x<=10000
    a, b, c, d =3.04218309, 1.17632271, 1.32270521, 0.17300204
    y = np.log(x)
    return np.exp(a * np.exp(-b * y) + c + d/y)   
         
        

def simul(fr,span,no,s2n=3.,seed=1949):
    '''
    t, v, er = simul(no,s2n=3.)
    Input:
    no - number of observations
    s2n - signal-to-noise {amplitude units]
    Output:
    t[no],v[no],er[no] simulated times, values end errors of a light curve
    '''
    rng = np.random.default_rng(seed)
    t=np.zeros(no)
    v=np.zeros(no)
    e=np.zeros(no)
    e[:] = np.exp(-rng.random(no))
    v[:] = np.random.normal(size=no)/s2n
    
    for k in np.arange(no):
       t[k] = span*np.sin(k*1.)**4
    t[-1] =span
    for k in np.arange(no):
       dph  = fr*t[k]
       omt  = np.pi*2*(dph-np.floor(dph))
       v[k]+= 11.90+0.434*np.cos(omt+4.72)+ \
           0.237*np.cos(2.*omt+0.741)
    return t,v,e
    
# a wrapper for the totals function

def test_per(nam):

    # get simulated data

    a = np.loadtxt(nam,comments="!")
    t, v, e = a[:,0],a[:,1],a[:,2]
    print("Close plots to see the next demo")
    # calculate the periodogram
    th,fr,frmax=amhw(t,v,e,2.8,0.0001,fr0=2.7,nh2=5) 
    # make a plot
    plper("amhw",frmax, t, v, fr, th)

    # calculate the periodogram
    th,fr,frmax=pspw(t,v,e,2.8,0.0001,fr0=2.7) 
    # make a plot
    plper("pspw",frmax, t, v, fr, th)

    # calculate the periodogram
    th,fr,frmax=atrw(t,v,e,2.8,0.0001,nh2=5,fr0=2.7) 
    print("This example of atrw is poor as data contain no sharp eclipses")
    # make a plot
    plper("atrw",frmax, t, v, fr, th)

    # calculate the periodogram
    th,fr,frmax=aovw(t,v,e,2.8,0.0001,nh2=5,fr0=2.7) 
    # make a plot
    plper("aovw",frmax, t, v, fr, th)
    
    # calculate the periodogram
    th,fr,frmax=f_mw(t,v,e,2.8,0.0001,fr0=2.7) 
    # make a plot
    plper("f_mw",frmax, t, v, fr, th)
    
    # calculate the periodogram
    th,fr,frmax=lomw(t,v,e,2.8,0.0001,fr0=2.7) 
    # make a plot
    plper("lomw",frmax, t, v, fr, th)
    


def aovdrv(method, tin, vin, win, npar, ncov, fr0, frs,nfr, sptype):
    """
    Interface routine to collect parameters from various aov methods and use the wrap method
    to pass the data to the fortran routine

    Parameters
    ----------
    method : string
        the algorithm (eg AOVW) to be used
    tin : numpy array of size (n*1)
        Times
    vin : numpy array of size (n*1)
        values 
    win : numpy array of size (n*1)
        errors
    nh2 : integer
        number of pphase bins
    ncov : integer
        number of coverages
    fr0 : float
        start frequency
    frs : float
        size of frequency steps
    nfr : integer
        number of frequency bins
    sptype : string
        spectrum type.

    Returns tuple
    -------
    th - numpy array (1*n) of powers 
    np0 - integer - number of parameters 
    stat - integer - status
    vr - float - variance of observations

    """

    _,stdout,stderr,parms,return_arrays=wrap("aovdrv",[method,npar,ncov,fr0,frs,nfr,sptype],[tin,vin,win])

    th=np.array([float(r) for r in return_arrays[0]])
    return [th,int(parms[0]),int(parms[1]),float(parms[2])]

def wrap(nam,parms,files):
    """
    This routine generates two temporary files - one for passing data to the fortran program
    and one to retrieve it. It then populates the first file using the three calling
    parameters and then invokes a subprocess.

    Parameters
    ----------
    nam : string 
        The name of the fortran subroutine to be invoked
    parms : list of parameters to be passed. The number and type will vary according to
             <name>
    files : list of  1-dimensional arrays.  
            The number and type will vary according to   <nam>
        

    Returns
    -------
    list consisting of:
        Status code - integer returned from the fortran routine
        wrap_info - string holding all the redirected output from stdout and stderr.
        files - :  list of datasets. Each dataset is a list of 1-dimensional arrays.  
                   The number and type will vary according to   <name>

    """
    import subprocess
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile('w+t',delete=False) as f:
        #print("Input filename is: ",f.name)
        with NamedTemporaryFile('w+t',delete=False) as g:
            #print("Output filename is: ",g.name)
            f.write (nam+"\n")
            f.write ('"'+g.name+'"\n')
            for v in parms:
                f.write(str(v)+",")
            f.write("\n")
            for fil in files:
                f.write(f"{str(len(fil))} \n")
                for i in range(len(fil)):
                    f.write(str(fil[i])+",")
                    f.write('\n')
            g.write("Fortran program did not complete\n")
    
            
    path = os.path.dirname(os.path.abspath(__file__))

    if path not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep+str(path)
    oldpath=os.environ["PATH"] 

    # As it is not possible to statically link Intel Fortran on windows
    # this code temporarily sets the path environment variable to 
    # how it was when aovwrapper was compiled
    if os.name=='nt':
        try:
            with open(os.sep.join([path,"path.dat"]),"r") as h:
                os.environ["PATH"]=h.read()[5:]+ os.pathsep+str(path)
        except:
            pass
    try:
        result=subprocess.run('aovwrapper '+f.name,capture_output=True,shell=True)
    except:
        os.environ["PATH"]=oldpath
        raise Exception("Failed to run aovwrapper")   
    os.environ["PATH"]    =oldpath  
    if result.returncode>0:

        print(f"Wrap subprocess.run error {result.returncode}\n\nstdout\n\n{result.stdout}\n\nstderr\n\n{result.stderr}\n")
        raise Exception("Wrap error")
    g2=open(g.name)
    reader=csv.reader(g2,delimiter=' ',skipinitialspace=True)
    nparameters=int(reader.__next__()[0])
    narrays=int(reader.__next__()[0])
    parms=[]
    for i in range(nparameters):
        parms+=[reader.__next__()[0]]
    return_arrays=[]
    for i in range(narrays):
        newarray=[]
        nitems=int(reader.__next__()[0])
        for n in range(nitems):
            newarray+=[reader.__next__()[0]]
        return_arrays+=[newarray]
    g2.close()    
    os.remove(f.name)
    os.remove(g.name)
    
    return [result.returncode,result.stdout,result.stderr,parms,return_arrays]
    
    
    

#==============================================  
# if you run this script from the command line this will
# run. It will not run if you import this file.
  
if __name__ == "__main__":
    

    test_per("test.dat")
    
