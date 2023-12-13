"""
This is a python GUI program to display and analyse a time series. It uses routines developed by 
Alex Schwarzenberg-Czerny to compute periodograms and pre-whitening (and alternatively a Fourier transform 
routine that is included).

It can either by called using a command line interface or as a function call:
    The function call (display_light_curve) is completely general and takes an astropy timeseries 
    object as its parameter (with specific metadata - see function definition)

    The command line is designed to process generic data and accepts 2 parameters (filename and path); 
    if these are not provided a dialog is initiated to select the values. The last path and format 
    used are saved in a .ini file located in the users home directory.

In addition there is a test capability that generates a synthetic light curve of two sine waves plus noise.
This is initiated by changing the value of the variable "test" to 1 and then
running the program. 

Prerequisites:
    This is python code using standard modules that
    are included in the Anaconda distribution. 

    It requires the implementation of the "Python GUI for F95 AOV periodogram routines"
    (https://users.camk.edu.pl/alex/) and calls the pyaov package which is included.
    This code has been designed to be portable and should work on
    Linux, Windows and Mac. It has been tested using Python 3.8 under
    Ubuntu 18.04 and also under Windows 10.

Coding notes
    This code is based on the QT (version QT5) widget toolkit
    (https://doc.qt.io/qtforpython/). This has been chosen rather than
    Tkinter because it is more portable and provides improved functionality
    (https://dev.to/amigosmaker/python-gui-pyqt-vs-tkinter-5hdd)
    A useful book is https://pythoncourses.gumroad.com/l/pysqtsamples
    
    Inevitably the structure of the program is driven by the GUI. display_light_curve
    initiates an instance of Window which will initiate an instance of Window1. Window1 
    can initiate an instance of Window2 which in turn can initiate an instance of Window3.
    
    
    



Change history

20/10/2022  Keith Inight        Initial version


"""
import sys
import os
from PyQt5.QtWidgets import QDialog, QApplication, QVBoxLayout, QHBoxLayout,\
    QPushButton, QFileDialog, QSpinBox, QLabel, QCheckBox, QButtonGroup, QRadioButton
from PyQt5.QtCore import QCoreApplication, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Chebyshev
from astropy.convolution import convolve, Box1DKernel
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel as W
from sklearn.exceptions import ConvergenceWarning
from scipy import optimize
from warnings import simplefilter
from numba import jit
import pprint
from ..Timeseries import pyaov
from bokeh.plotting import figure

pp = pprint.PrettyPrinter(width=130)

@jit(nopython=True)
def fourier_transform(nout, fmin, fmax, times, mags):
    """      
     Purpose
         To compute the Fourier transform (power curve) of a set of data. 

         Based on algorithm  - https://ui.adsabs.harvard.edu/abs/1975Ap%26SS..36..137D/abstract
     Parameters 

        nout:    int                   Number of frequency bins in output
        fmin :   float                 The minimum frequency of the frequency bins ( the frequency
                                           unit is the inverse of the unit of time)       
        fmax :   float                 The maximum frequency of the frequency bins ( the frequency
                                           unit is the inverse of the unit of time)
        times:   np.ndarray of float   values of time variable
        mags:    np.ndarray of float   values of magnitude variable

     Returns 
        A list of 4 objects:-
        freqs   np.ndarray of float    value of frequency for each bin
        periods np.ndarray of float    value of period for each bin (in same unit as times)
        FF      np.ndarray of float    amplitude squared for each bin
        GG      np.ndarray of powers   normalized amplitude squared for each bin
     History
         11/2/20 Keith Inight. Initial version
    """

    ndata = len(times)
    fr = np.zeros(nout)
    fi = np.zeros(nout)
    D = np.zeros(nout)
    G = np.zeros(nout)
    FF = np.zeros(nout)
    GG = np.zeros(nout)
    freqs = np.linspace(fmin, fmax, nout)
    periods = 1/freqs
    angular_freqs = (2*np.pi*freqs)
    for k in range(nout):
        for i in range(ndata):
            A = angular_freqs[k]*times[i]
            c = np.cos(A)
            s = np.sin(A)
            fr[k] += mags[i]*c
            fi[k] += mags[i]*s
            D[k] += c
            G[k] += s
        FF[k] = (fr[k]*fr[k]+fi[k]*fi[k])/(ndata**2)
        GG[k] = (D[k]*D[k]+G[k]*G[k])/(ndata**2)
    return [freqs, periods, FF, GG]


class Window3(QDialog):
    def __init__(self, parent=None, ts_detrended=[]):

        super(Window3, self).__init__(parent)
        self.ts_detrended = ts_detrended
        self.fpeak = self.ts_detrended.meta['fpeak']
        self.fpeak_err = self.ts_detrended.meta['fpeak_err']
        self.setWindowTitle("{} ".format(
            self.ts_detrended.meta['object_name']))
        # a figure instance to plot on
        self.figure = plt.figure(dpi=80)
        self.setGeometry(200, 200, 1000, 600)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__

        self.canvas3 = FigureCanvas(self.figure)

        self.toolbar3 = NavigationToolbar(self.canvas3, self)

        period_fold = 1/self.fpeak

        self.phase = self.ts_detrended.time.value/period_fold - \
            (self.ts_detrended.time.value/period_fold).astype('int')
        if len(self.ts_detrended['y3']) > 50:
            bins = 30
        else:
            bins = 10
        np.linspace(0, 1-1/bins, num=bins)
        self.bin_phase = []
        self.bin_y = []
        self.bin_y_err = []
        for bphase in np.linspace(0, 1-1/bins, num=bins):
            mask = (self.phase < bphase+1/bins) & (self.phase >= bphase)
            if len(self.ts_detrended['y3'][mask]) > 0:
                self.bin_phase += [bphase+0.5/bins]
                baverage, norm = np.average(self.ts_detrended['y3'][mask].value, weights=(
                    1/self.ts_detrended['y_err'][mask].value)**2, returned=True)
                self.bin_y += [baverage]
                self.bin_y_err += [1/np.sqrt(norm)]
        self.plot()

        # set the layout
        layout3 = QVBoxLayout()
        layout3.addWidget(self.toolbar3)
        layout3.addWidget(self.canvas3)
        self.setLayout(layout3)

    def closeEvent(self, event):
        plt.close()
        super(Window3, self).closeEvent(event)

    def plot(self):
        self.figure.clear()
        gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1])
        y3 = self.ts_detrended['y3'].value
        err = self.ts_detrended['y_err'].value

        # create an axis
        ax = plt.subplot(gs[0, 0])
        # plot data
        ax.errorbar(np.concatenate((self.phase, self.phase+1)),
                    np.concatenate((y3, y3)), yerr=np.concatenate((err, err)), fmt='.', ecolor='grey', elinewidth=.3)

        ax.set_title('f = {:.6f}'.format(self.fpeak)+' d$^{-1}$' +
                     '   P={:.6f} $\pm$ {:.8f}'.format(24/self.fpeak, 24*self.fpeak_err/((self.fpeak)**2)) +
                     ' h ({:.6f}  $\pm$ {:.10f} days)'.format(1/self.fpeak, self.fpeak_err/(self.fpeak**2)))

        ax.set_xlabel('Phase')
        ax.set_ylabel(str(self.ts_detrended['y3'].unit))
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymax, ymin)
        # refresh canvas
        ax1 = plt.subplot(gs[0, 1])
        bp = np.array(self.bin_phase)
        ax1.errorbar(np.concatenate([bp, bp+1]),
                     self.bin_y + self.bin_y,
                     self.bin_y_err + self.bin_y_err, fmt='.', ecolor='grey', elinewidth=.3)

        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim(ymax, ymin)
        ax1.set_title('Binned plot')
        ax1.set_xlabel('Phase')
        ax2 = plt.subplot(gs[1, 0:2])
        ax2.text(0, 0.95, pp.pformat({i: self.ts_detrended.meta[i] for i in self.ts_detrended.meta if i in
                                      ['frange', 'crop', 'prewhiten', 'trend', 'Gaussian_Process_regression', 'method','peak_type', 'nh2', 'fpeak']}), va='top')
        ax2.axis('off')

        self.canvas3.draw()


class Window2(QDialog):
    def __init__(self, parent=None, ts_detrended=None):

        super(Window2, self).__init__(parent)

        self.ts_detrended = ts_detrended
        self.setWindowTitle("{} ".format(ts_detrended.meta['object_name']))
        # a figure instance to plot on
        self.figure = plt.figure(dpi=80)
        self.setGeometry(150, 150, 800, 400)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas1 = FigureCanvas(self.figure)
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        self.times = self.ts_detrended.time.value-self.ts_detrended.time.value.min()
        self.y = self.ts_detrended['y3'].value
        self.y_err = self.ts_detrended['y_err'].value
        method = self.ts_detrended.meta['method']
        self.nh2 = self.ts_detrended.meta['nh2']
        self.binsize = self.ts_detrended.meta['binsize']

        fmin, fmax, fbins = self.ts_detrended.meta['frange'][-1]


        if method == "amhw":    
            self.per, self.freq, self.fpeak = pyaov.amhw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin, nh2=self.nh2)

        elif method == "pspw":

            self.per, self.freq, self.fpeak = pyaov.pspw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin)

        elif method == "atrw":

            self.per, self.freq, self.fpeak = pyaov.atrw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin, nh2=self.nh2)

        elif method == "aovw":

            self.per, self.freq, self.fpeak = pyaov.aovw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin, nh2=self.nh2)

        elif method == "f_mw":

            self.per, self.freq, self.fpeak = pyaov.f_mw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin)

        elif method == "lomw":

            self.per, self.freq, self.fpeak = pyaov.lomw(
                self.times, self.y, self.y_err, fmax, (fmax-fmin)/fbins, fr0=fmin)

        else:
            print("Invalid method ", method)
            raise ValueError("Invalid method "+str(method))

        self.plot(self.freq, self.per, title=method)

        self.Fold_button = QPushButton("Fold")
        self.Fold_button.clicked.connect(self.Phase_fold)

        self.Fold_half_button = QPushButton("Fold "+"\u00BD"+"f")
        self.Fold_half_button.clicked.connect(self.Phase_fold_half)

        self.Recalculate_button = QPushButton("Recalculate")
        self.Recalculate_button.clicked.connect(self.Recalculate)

        self.Prewhiten_button = QPushButton("Prewhiten")
        self.Prewhiten_button.clicked.connect(self.Prewhiten)

        self.row_layout = QHBoxLayout()
        self.row_layout.addWidget(self.toolbar1)
        self.row_layout.addWidget(self.Fold_button)
        self.row_layout.addWidget(self.Fold_half_button)
        self.row_layout.addWidget(self.Recalculate_button)
        self.row_layout.addWidget(self.Prewhiten_button)
        # set the layout
        layout1 = QVBoxLayout()
        layout1.addLayout(self.row_layout)
        layout1.addWidget(self.canvas1)
        self.setLayout(layout1)

    def closeEvent(self, event):

        plt.close()
        super(Window2, self).closeEvent(event)

    def Phase_fold(self):
        def sine_fit(x, a, b, c):
            return a * np.sin(b * x+c)
        axes = self.canvas1.figure.axes[0]
        freq_lo, freq_hi = axes.get_xlim()

        mask = (self.freq > freq_lo) & (self.freq < freq_hi)
        if (self.fpeak > freq_lo) & (self.fpeak< freq_hi):
                self.ts_detrended.meta['fpeak']=self.fpeak
        else:
                self.ts_detrended.meta['fpeak'] = self.freq[mask][self.per[mask].argmax()]  # seed value
        if self.ts_detrended.meta['peak_type']=='Parabola': 
            xm, fm, dx = pyaov.peak(self.per[mask])
    
            self.ts_detrended.meta['fpeak'] = self.freq[mask][0]+xm*self.binsize
            self.ts_detrended.meta['fpeak_err'] = dx*self.binsize
        else:
            print("OTHER")
            
            try:
                print("fpeak",self.ts_detrended.meta['fpeak'])
                self.params, self.params_covariance = optimize.curve_fit(sine_fit, self.ts_detrended.time.value, self.ts_detrended['y3'].value,
                                                                         p0=[1, 2*np.pi*self.ts_detrended.meta['fpeak'], 0])
                print(self.params)
                print(self.params_covariance)
            except:
                print('failed')
                self.params = np.array([0, 2*np.pi*self.ts_detrended.meta['fpeak'], 0])
                self.params_covariance = np.array(
                    [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
            
  
            self.ts_detrended.meta['fpeak'] = self.params[1]/(2*np.pi)
            self.ts_detrended.meta['fpeak_err'] =     np.sqrt(self.params_covariance[1,1])/(2*np.pi)
            print("output",self.ts_detrended.meta['fpeak'])
        
        main3 = Window3(self, ts_detrended=self.ts_detrended)
        main3.show()

    def Phase_fold_half(self):

        axes = self.canvas1.figure.axes[0]
        freq_lo, freq_hi = axes.get_xlim()
        mask = (self.freq > freq_lo) & (self.freq < freq_hi)
        xm, fm, dx = pyaov.peak(self.per[mask])
        self.ts_detrended.meta['fpeak'] = (
            self.freq[mask][0]+xm*self.binsize)*0.5
        self.ts_detrended.meta['fpeak_err'] = dx*self.binsize
        main3 = Window3(self, ts_detrended=self.ts_detrended)
        main3.show()

    def Recalculate(self):

        axes = self.canvas1.figure.axes[0]
        freq_lo, freq_hi = axes.get_xlim()
        fbins = self.ts_detrended.meta['frange'][-1][2]

        self.ts_detrended.meta['frange'] += [[freq_lo, freq_hi, fbins]]
        
        self.ts_detrended.meta['binsize'] = (freq_hi -
                              freq_lo)/fbins
        
        main2a = Window2(self, ts_detrended=self.ts_detrended)
        main2a.show()

    def Prewhiten(self):
        axes = self.canvas1.figure.axes[0]
        freq_lo, freq_hi = axes.get_xlim()
        per_lo, per_hi = axes.get_ylim()
        mask = (self.freq > freq_lo) & (self.freq < freq_hi) & (
            self.per > per_lo) & (self.per < per_hi)

        xm, fm, dx = pyaov.peak(self.per[mask])

        foundfreq = self.freq[mask][0]+xm*self.binsize
        # already got this frequency in the list

        frout,dfrout,valout=pyaov.prew(self.times, self.y, self.y_err, foundfreq, nh2=11)
        print('frout={}  dfrout={}'.format(frout,dfrout))

        if frout in self.ts_detrended.meta['prewhiten']:
            return
        
        self.ts_detrended.meta['fpeak'] =frout
        
        self.ts_prewhite = self.ts_detrended
             
        self.ts_prewhite.meta['prewhiten'] += [foundfreq]
        print(self.ts_prewhite.meta['prewhiten'])

        self.ts_prewhite['y3'] =  valout*self.ts_detrended['y3'].unit 

        self.ts_prewhite.meta['frange'] = [self.ts_detrended.meta['frange'][0]]
        main2b = Window2(self, ts_detrended=self.ts_prewhite)
        main2b.show()

    def plot(self, freq, per, title=''):

        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)
        # plot data
        ax.plot(freq, per, '-', linewidth=0.5)
        ax.set_xlabel('Frequency(1/D)')
        ax.set_ylabel('Power')
        ax.set_title(title)
        # refresh canvas
        self.canvas1.draw()


class Window1(QDialog):
    def __init__(self, parent=None, ts_cropped=None):
        super(Window1, self).__init__(parent)
        self.ts_cropped = ts_cropped
        self.setWindowTitle("{} ".format(
            str(self.ts_cropped.meta['object_name'])))
        # a figure instance to plot on
        # self.ts_cropped.sort('time')

        self.figure = plt.figure(dpi=80)
        self.setGeometry(100, 100, 1000, 400)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent

        #self.toolbar = NavigationToolbar(self.canvas, self)

        self.label1 = QLabel("Trend")
        self.label1.alignment()
        self.label1.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.trendbox = QSpinBox()
        self.trendbox.setGeometry(0, 0, 30, 20)
        self.trendbox.setValue(0)
        self.trendbox.setMinimum(0)
        self.trendbox.valueChanged.connect(self.redraw)
        # -------------------------
        self.label2 = QLabel("Boxcar")
        self.label2.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.label2.setToolTip("This sets the length of the boxcar")

        self.boxcarbox = QSpinBox()
        self.boxcarbox.setGeometry(0, 0, 30, 20)
        self.boxcarbox.setValue(1)
        self.boxcarbox.setMinimum(1)
        self.boxcarbox.setSingleStep(2)
        self.boxcarbox.valueChanged.connect(self.redraw)

        # -------------------------
        self.Gaussian_Process_Regression = QCheckBox(
            "Gaussian Process Regression")
        self.Gaussian_Process_Regression.setLayoutDirection(Qt.RightToLeft)
        self.Gaussian_Process_Regression.stateChanged.connect(self.redraw)
        # -------------------------

        self.label3 = QLabel("nh")
        self.label3.alignment()
        self.label3.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.nhbox = QSpinBox()
        self.nhbox.setGeometry(0, 0, 30, 20)
        self.nhbox.setValue(2)
        self.nhbox.setMinimum(2)

        # -------------------------
        self.row_layout = QHBoxLayout()
        self.row_layout.addWidget(self.label1)
        self.row_layout.addWidget(self.trendbox)
        self.row_layout.addWidget(self.label2)
        self.row_layout.addWidget(self.boxcarbox)
        self.row_layout.addWidget(self.Gaussian_Process_Regression)
        self.row_layout.addWidget(self.label3)
        self.row_layout.addWidget(self.nhbox)

        self.fpeaktype=QButtonGroup()
        self.fpeaktype.setExclusive(True)

        self.radio = QRadioButton("Parabola")
        self.fpeaktype.addButton(self.radio, 0)
        self.radio.setChecked(True)
        self.radio.setLayoutDirection(Qt.RightToLeft)
        self.row_layout.addWidget(self.radio)
        self.row_layout.setAlignment(self.radio, Qt.AlignVCenter)
        self.radio = QRadioButton("Sine")
        self.fpeaktype.addButton(self.radio, 1)
        self.radio.setChecked(False)
        self.radio.setLayoutDirection(Qt.RightToLeft)
        self.row_layout.addWidget(self.radio)
        self.row_layout.setAlignment(self.radio, Qt.AlignVCenter)
        
        #self.fpeaktype.buttonClicked.connect(self.on_draw_per)



        #  Define the second row of widgets

       # self.nhbox.valueChanged.connect(self.nh2)
       # self.nh2=2

        self.listOfChoices = ["amhw", "pspw", "atrw",
                              "aovw", "f_mw", "lomw"]
        self.choice_tooltips = {"amhw": "Returns multiharmonic AOV periodogram,\n"
                                        "obtained by fitting data with a series of trigonometric polynomials.\n"
                                        "For default nh2=3 this is Lomb-Scargle periodogram\n"
                                        "corrected for constant shift",
                                "pspw": "Returns Power Spectrum. Principal use for calculation of window function.\n"
                                        "For periodogram we recommend amhv routine with nh2=2 instead",
                                "atrw": "Returns AOV periodogram for search of transits/eclipses.\n"
                                        "For this purpose data are fitted with a box-like \n"
                                        "(inverted top-hat) function",
                                "aovw": "Analysis of variance taking account of weights\n"
                                        "i.e. observastional errors",
                                "f_mw": "Returns generalized Lomb-Scargle periodogram\n"
                                        "of data, using the original Ferraz-Mello (1981) algorithm.",
                                "lomw": "Returns AOV periodogram, obtained by phase-folding and binning of data."}
        self.row2_layout = QHBoxLayout()

        self.radioGroup = QButtonGroup()
        self.radioGroup.setExclusive(True)

        for i, row in enumerate(self.listOfChoices):
            self.radio = QRadioButton(row)
            self.radio.setToolTip(self.choice_tooltips[row])
            self.radioGroup.addButton(self.radio, i)
            # if i==0:
            #   self.radio.setChecked(True)

            self.row2_layout.addWidget(self.radio)
            self.row2_layout.setAlignment(self.radio, Qt.AlignVCenter)
        self.radioGroup.buttonClicked.connect(self.on_draw_per)

        # set the top level layout
        layout = QVBoxLayout()
        layout.addLayout(self.row_layout)
        layout.addLayout(self.row2_layout)
        layout.addWidget(self.canvas)

        self.setLayout(layout)
        self.plot()
        simplefilter('ignore', category=ConvergenceWarning)

    def on_draw_per(self, meth):
        print(self.fpeaktype.checkedButton().text())
        self.ts_cropped.meta['peak_type']=self.fpeaktype.checkedButton().text()
        self.ts_cropped.meta['method'] = meth.text()
        main2 = Window2(self, ts_detrended=self.get_detrended_subset())
        main2.show()

    def closeEvent(self, event):

        plt.close()
        super(Window1, self).closeEvent(event)

    def get_detrended_subset(self):
        # calculate best fit trend line
        global dummy
        self.ts_detrended = self.ts_cropped.copy()

        self.ts_cropped.meta['trend'] = self.trendbox.value()
        self.ts_cropped.meta['boxcar'] = self.boxcarbox.value()
        self.ts_cropped.meta['nh2'] = self.nhbox.value()
        self.ts_cropped.meta['Gaussian_Process_Regression'] = self.Gaussian_Process_Regression.isChecked(
        )

        c = Chebyshev.fit(self.ts_detrended.time.value,
                          self.ts_detrended['y'].value, deg=self.trendbox.value())
        self.ts_detrended['y1'] = self.ts_detrended['y'] - \
            c(self.ts_detrended.time.value)*self.ts_detrended['y'].unit

        if self.boxcarbox.value() > 1:
            box_kernel = Box1DKernel(self.boxcarbox.value())
            smoothed = convolve(self.ts_detrended['y1'].value, box_kernel)
            self.ts_detrended['y2'] = self.ts_detrended['y1'] - \
                smoothed*self.ts_detrended['y1'].unit
            for i in range(int(self.boxcarbox.value()/2)):
                self.ts_detrended['y2'][i] = 0
                self.ts_detrended['y2'][len(self.ts_detrended['y2'])-1-i] = 0
        else:
            self.ts_detrended['y2'] = self.ts_detrended['y1']

        # self.ts_detrended.sort('time')

        if self.Gaussian_Process_Regression.isChecked():
            kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1, 1e2)) + \
                W(noise_level=1e-3, noise_level_bounds=(1e-06, 1))
            gp = GaussianProcessRegressor(
                kernel=kernel, alpha=self.ts_cropped['y_err'].value ** 2, n_restarts_optimizer=3)
            time2d = np.atleast_2d(self.ts_detrended.time.value).T
            gp.fit(time2d, self.ts_detrended['y2'].value)
            mag_pred, sigma_mag_pred = gp.predict(time2d, return_std=True)
            self.ts_detrended['y3'] = self.ts_detrended['y2'] - \
                mag_pred * self.ts_detrended['y2'].unit
        else:
            self.ts_detrended['y3'] = self.ts_detrended['y2']
        return self.ts_detrended

    def Fourier(self):
        self.ts_cropped.meta['method'] = 'Fourier'
        main2 = Window2(self, ts_detrended=self.get_detrended_subset())
        main2.show()

    def redraw(self):
        # Update display to include the updated trend/boxcar parameters to be
        # displayed
        self.plot()

    def plot(self):

        self.ts_detrended = self.get_detrended_subset()

        xdata = self.ts_cropped.time.value
        # calculate the "trend" that has been used
        ydata = self.ts_detrended['y']-self.ts_detrended['y3']

        self.figure.clear()
        # create an axis
        ax = self.figure.add_subplot(111)
        # plot data
        ax.errorbar(self.ts_cropped.time.value,
                    self.ts_cropped['y'].value, yerr=self.ts_cropped['y_err'].value, fmt='.', color='magenta')
        ax.plot(xdata, ydata, '-k')
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymax, ymin)

        ax.set_xlabel('MJD')
        ax.set_ylabel('Magnitude')
        ax.set_title(str(self.ts_cropped.meta['object_name']))
        # refresh canvas

        self.canvas.draw()

class Window(QDialog):
    def __init__(self, ts, parent=None):
        super(Window, self).__init__(parent)
        self.ts = ts
        self.setWindowTitle("{} ".format(str(self.ts.meta['object_name'])))
        # a figure instance to plot on

        self.figure = plt.figure(dpi=80)
        self.setGeometry(50, 50, 1000, 400)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.analyse_button = QPushButton("Analyse")
        self.analyse_button.clicked.connect(self.Analyse)

        self.row_layout = QHBoxLayout()
        self.row_layout.addWidget(self.toolbar)
        self.row_layout.addWidget(self.analyse_button)

        # set the top level layout
        layout = QVBoxLayout()
        layout.addLayout(self.row_layout)
        layout.addWidget(self.canvas)

        self.setLayout(layout)
        self.plot()

    def closeEvent(self, event):

        plt.close()
        super(Window, self).closeEvent(event)

    def get_subset(self):

        try:
            axes = self.canvas.figure.axes[0]
            mjd1, mjd2 = axes.get_xlim()
            y_hi, y_lo = axes.get_ylim()
            mask = (self.ts.time.value > mjd1) &\
                (self.ts.time.value < mjd2) &\
                (self.ts['y'].value > y_lo) &\
                (self.ts['y'].value < y_hi)
            self.new_ts = self.ts.copy()[mask]
            self.new_ts.meta['crop'] = [mjd1, mjd2, y_lo, y_hi]
        except:
            self.new_ts = self.ts.copy()
        return self.new_ts

    def Analyse(self):

        main1 = Window1(self, self.get_subset())
        main1.show()

    def redraw(self):
        # Update display to include the updated trend/boxcar parameters to be
        # displayed
        self.plot()

    def plot(self):
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)
        # plot data
        ax.errorbar(self.ts.time.value, self.ts['y'].value,
                    yerr=self.ts['y_err'].value, fmt='.', color='magenta')
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymax, ymin)
        ax.set_xlabel(str(self.ts.time.format))
        ax.set_ylabel(str(self.ts['y'].unit))
        ax.set_title(str(self.ts.meta['object_name']))
        # refresh canvas
        self.canvas.draw()


def display_light_curve(original_ts):
    """
    Parameters
    ----------
    original_ts : Timeseries object
        This is a standard astropy time series object with (dimensioned) parameters
        time,y and y_err. In addition it includes a metadata dictionary (meta) with
        parameters
            'object path': path to folder where time series data is held
            'object_name': name of object
            'frange': list of lists. The lowest level list has three values applicable
                to periodograms ; lowest frequency, highest frequency and number of bins
           
    Returns
    -------
    None
    """
    ts = original_ts.copy()
    ts.sort('time')
    if 'object_name' not in ts.meta:
        ts.meta['object_name'] = ''
    ts.meta['crop'] = []
    if 'prewhiten' not in ts.meta:
        ts.meta['prewhiten'] = []
    if 'frange' not in ts.meta:
        period_low = 0.5  # shortest period evaluated (hours)
        period_hi = 30  # longest period evaluated
        fmin = 24/period_hi
        fmax = 24/period_low
        fbins = 150000
        ts.meta['frange'] = [[fmin, fmax, fbins]]
    ts.meta['binsize'] = (ts.meta['frange'][0][1] -
                          ts.meta['frange'][0][0])/ts.meta['frange'][0][2]
    
    # Cannot change order of these
    app = QApplication(sys.argv)
    main = Window(ts)
    main.show()
    app.exec_()

def getanalysis(data):
    from astropy.time import Time
    from astropy.timeseries import TimeSeries
    from astropy.units import cds
    from astropy import units as u
    from astropy.io import ascii, fits  

    ts_filename='lightcurve_data.csv'
    data.to_csv(ts_filename,index=False)

    ts_path = os.path.join(os.getcwd(),ts_filename)

    flux_unit=u.def_unit("Flux")
    no_unit=u.def_unit("Unspecified")
    
    '''
    _, appname = os.path.split(sys.argv[0])

    default_file = os.sep.join([str(Path.home()), 'KI app parameters', appname[:-3]+'.def'])
    '''
    
    data.sort_values('hjd',inplace=True)
    xd=data['hjd']-2400000.5

    yd=data['mag']
    ed=data['magerr']

    y_unit=u.mag

    filters = list(set(data['filtercode']))
    first_filter = filters[0]
    med_first_filter = data.query('filtercode==@first_filter')['mag'].median()
    
    for filt in filters:
        med_filtered = data.query('filtercode==@filt')['mag'].median()
        data.loc[data['filtercode'] == filt,'mag'] += (med_first_filter-med_filtered)
    yd=data['mag']

    ts_time = Time(xd.to_numpy(), format='mjd')
    ts_timeseries = TimeSeries(time=ts_time, meta={'object path': ts_path, 'object_name': ts_filename, 'frange': [[2, 60, 1500000]]})
    
    ts_timeseries.add_column(yd.to_numpy()*y_unit, name='y')
    ts_timeseries.add_column(ed.to_numpy()*cds.mag, name='y_err')

    display_light_curve(ts_timeseries)

    '''
    f = open(default_file, 'w')
    f.write(ts_path)
    f.close()
    '''
    
def getpowerspectrum(data):
    from astropy.time import Time
    from astropy.timeseries import TimeSeries
    from astropy.units import cds
    from astropy import units as u
    
    file_name='lightcurve_data.csv'
    data.to_csv(file_name,index=False)
    
    file_path=os.path.join(os.getcwd(),file_name)
    
    data.sort_values('hjd',inplace=True)
    xd=data['hjd']-2400000.5
    yd=data['mag']
    ed=data['magerr']
    
    y_unit=u.mag
    
    filters = list(set(data['filtercode']))
    first_filter = filters[0]
    med_first_filter = data.query('filtercode==@first_filter')['mag'].median()
    
    for filt in filters:
        med_filtered = data.query('filtercode==@filt')['mag'].median()
        data.loc[data['filtercode'] == filt,'mag'] += (med_first_filter-med_filtered)
    yd=data['mag']
    
    start_freq=0
    final_freq=50
    step_size=(final_freq-start_freq)/250000

    th,freqs,frmax=pyaov.amhw(xd,yd,ed,final_freq,step_size)
    
    plot=figure(width=400,height=400,title='ZTF Power Spectrum',x_axis_label=r'\[\text{Frequency [days}^{-1}]\]',y_axis_label=r'\[\text{Power}\]')
    plot.line(x=freqs,y=th,color='black',line_width=0.5)
    
    return plot