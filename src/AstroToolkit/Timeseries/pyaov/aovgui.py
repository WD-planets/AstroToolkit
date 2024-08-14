#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
GUI for AOV Period Search in Light Curves

Author: Alex Schwarzenberg-Czerny, alex@camk.edu.pl

To be done:  
   -add flattening of transit light curves by ort fit with
    clipping of outstanding points.
   -use showFullScreen() on your widget
   -define ReDo button
   -define starlist input
"""

# Packages
import sys
import os
import numpy as np
#import copy
#from scipy import *

import astropy.io.ascii as asc
from astropy.io import ascii, fits
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon

from pylab import getp

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.backend_bases import MouseButton

from PyQt5.QtWidgets import QApplication, QVBoxLayout, QHBoxLayout,\
    QPushButton, QFileDialog, QLabel, QCheckBox, QMainWindow, QAction,\
    QWidget, QDesktopWidget, QLineEdit, QButtonGroup, QRadioButton, QSlider, QMessageBox

#import pyaovv2 as pyaov
import pyaov


class My_Lc:  # Light curve data
    def __init__(self):  # Initialize program
        self.fname = ""
        self.xd = np.array([])
        self.yd = np.array([])
        self.ed = np.array([])
        self.nact = 0
        self.sw = 0.
        self.av = 0.
        self.vr = 1.
        self.epoch = 0.

    def __str__(self):
        s = ''+self.fname
        n = self.xd.size
        print('n=', n)
        if n < 5:
            l = np.arange(n)
        else:
            l = np.array([0, 1, -2, -1])
        for i in l:
            if n >= 5 and i == -2:
                s += ' ...'
            s += ' '+str(self.xd[i])
        for i in l:
            if n >= 5 and i == -2:
                s += ' ...'
            s += ' '+str(self.yd[i])
        for i in l:
            if n >= 5 and i == -2:
                s += ' ...'
            s += ' '+str(self.ed[i])
        s += ' '+str(self.nact)
        s += ' '+str(self.sw)
        s += ' '+str(self.av)
        s += ' '+str(self.vr)
        s += ' '+str(self.epoch)
        return s


def Cp_Lc(a, b):
    b.fname = str(a.fname)
    b.xd = np.array(a.xd)
    b.yd = np.array(a.yd)
    b.ed = np.array(a.ed)
    b.nact = (a.nact)
    b.sw = (a.sw)
    b.av = (a.av)
    b.vr = (a.vr)
    b.epoch = (a.epoch)


class My_Frm:  # Basic parameters of a plot frame
    def __init__(self):  # Initialize program
        self.xd = np.array([])
        self.yd = np.array([])
        self.xr = np.zeros(3, dtype=float)  # xlow,xup,xstep
        self.nx = 0
        self.yr = np.zeros(2, dtype=float)  # ylow,yup
        self.replt = False                # replot figure
        self.pxy = np.zeros(2, dtype=float)  # point coords
        self.pplt = False                 # plot point?
        self.ax = None                    # ax handle(subplot)
        self.plot = None                  # plot routine

    def __str__(self):
        s = ''
        n = self.xd.size
        if n < 5:
            l = np.arange(n)
        else:
            l = np.array([0, 1, -2, -1])
        for i in l:
            if n >= 5 and i == -2:
                s += ' ...'
            s += ' '+str(self.xd[i])
        for i in l:
            if n >= 5 and i == -2:
                s += ' ...'
            s += ' '+str(self.yd[i])
        for i in self.xr:
            s += ' '+str(i)
        s += ' '+str(self.nx)
        for i in self.yr:
            s += ' '+str(i)
        s += ' '+str(self.replt)
        for i in self.pxy:
            s += ' '+str(i)
        s += ' '+str(self.pplt)
        s += ' '+str(self.ax)
        s += ' '+str(self.plot)
        return s


def Cp_Frm(a, b):
    b.xd = np.array(a.xd)
    b.yd = np.array(a.yd)
    b.xr = np.array(a.xr)
    b.nx = (a.nx)
    b.yr = np.array(a.yr)
    b.replt = (a.replt)
    b.pxy = np.array(a.pxy)
    b.pplt = (a.pplt)
    b.ax = a.ax
    b.plot = a.plot


class AppForm(QMainWindow):

    def __init__(self, parent=None):  # Initialize program
        QMainWindow.__init__(self, parent)
#        screen = QDesktopWidget().screenGeometry()
#        size = self.geometry()
#        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)
#        self.setGeometry(int(screen.width()*0.7), int(screen.height()*0.7), \
#               screen.width()/2, screen.height()/2)
        self.strold = ""
        self.marg = 0.07
        self.setWindowTitle('AOV Period Search')
        self.per0 = My_Frm()
        self.per0.plot = self.on_draw_per
        self.lc0 = My_Frm()
        self.lc0.plot = self.on_draw_lc

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.open_dialog()  # open file & set default values

    def open_dialog(self):  # Action for File>Open
        fd = QFileDialog(self)
        plik = str(fd.getOpenFileName()[0])
        fd.close()
        # input light curve
        self.lcd0 = self.input_lc(plik)
        self.lcd = My_Lc()
        Cp_Lc(self.lcd0, self.lcd)
        self.lcd.yd, self.lcd.nact, self.lcd.sw, self.lcd.av, self.lcd.vr = \
            pyaov.normalize(self.lcd.yd, self.lcd.ed, var=0.)

        # initialize frame parameters
        print("Time range ", self.lcd.xd[0], self.lcd.xd[-1])
        self.per0.xr[1], self.per0.xr[2], self.per0.xr[0] = \
            pyaov.fgrid(self.lcd.xd)
        self.per0.xr[1] = 1.5*self.per0.xr[1]
        # self.per0.xr[1] = 5*self.per0.xr[1] # KPIaddition
        self.per = My_Frm()
        Cp_Frm(self.per0, self.per)
        self.last_click = None
        self.nh2 = 3
        self.nh2box.setText("%3d" % self.nh2)

        self.lc0.xr[0], self.lc0.xr[1] = (0., 1.2)
        self.lc0.xr[2] = (self.lc0.xr[1]-self.lc0.xr[0])/1200
        self.lc = My_Frm()
        Cp_Frm(self.lc0, self.lc)

        # get periodogram
        self.on_draw_per(self.radioGroup.checkedId())
        # print(unicode(self.textbox.text()))

    def freeze_dialog(self):  # Action for File>Freeze & Open
        fd = QFileDialog(self)
        plik = str(fd.getOpenFileName()[0])
        fd.close()

        self.lcd0 = self.input_lc(plik)  # input light curve
        self.lcd = My_Frm()
        Cp_Frm(self.lcd0, self.lcd)
        self.lcd.yd, self.lcd.nact, self.lcd.sw, self.lcd.av, self.lcd.vr = \
            pyaov.normalize(self.lcd.yd, self.lcd.ed, var=0.)

        # freeze & get periodogram
        self.per0.xr = np.array(self.per.xr)
        self.lc0.xr = np.array(self.lc.xr)
        self.last_click = None
        self.on_draw_per(self.radioGroup.checkedId())
#        print(unicode(self.textbox.text()))

    def save_plot(self):  # Action for File>Save Plot
        file_choices = "PNG (*.png)|*.png"

       # path = unicode(QFileDialog.getSaveFileName(self, \
        #          'Save file', '', file_choices))
        path = QFileDialog.getSaveFileName(self,
                                           'Save file', '', file_choices)[0]

        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)

    def save_data(self):  # Action for File>Save Data
        file_choices = "DAT (*.dat)|*.dat"

       # path = unicode(QFileDialog.getSaveFileName(self, \
        #          'Save file', '', file_choices))
        path = QFileDialog.getSaveFileName(self,
                                           'Save file', '', file_choices)[0]
        if path:
            asc.write({'x': self.xd, 'y': self.yd, 'z': self.ed},
                      path, Writer=asc.NoHeader)
            self.statusBar().showMessage('Saved to %s' % path, 2000)

    def save_period(self):  # Action for File>Save Periodogram
        file_choices = "DAT (*.per)|*.per"

       # path = unicode(QFileDialog.getSaveFileName(self, \
        #          'Save file', '', file_choices))
        path = QFileDialog.getSaveFileName(self,
                                           'Save file', '', file_choices)[0]
        if path:
            asc.write({'x': self.per.xd, 'y': self.per.yd},
                      path, Writer=asc.NoHeader)
            self.statusBar().showMessage('Saved to %s' % path, 2000)

    def on_about(self):  # Action on About>Help
        msg = """ AOV Period Search of Light Curve
   To start this program e.g. type: python aovgui.py | tee aov.log. In this way messages would be written both to the screen and to the aov.log file. You may 
comment your log with NOTEs in the bottom box.
   * Use File>OPEN to input data in 2 or 3 column ascii (time,value[,error]);
# commented headers are permitted.
   * Top and bottom plots present respectively periodogram and data  folded with
the peak frequency (marked with dot)
   * Try pointing other peaks with mouse (within 0.1-0.9 y-scale) or type 
frequency to change data folding
   * Zoom-in/out the plots by clicking bottom/top frame on two locations to select the range. Zoom-in would cover the selected range while zoom-out would expand the x-scale asymmetricaly by the ratio of original and selected ranges...
   * ...or refine your frequency by fitting trig polynomial series. To obtain 
a realistic estimate of the frequency error, use REFINE when displaying moderate 
range of amh|aov|atr periodogram, centered on the line: median of the displayed  spectrum is used to estimate noise level.
   * individual observations may be rejected/re-included by clicking them on the 
LC plot. Note, this plot can also be zoomed. Enter 0. frequency to plot LC against time [time unit=data span]
   * PREWHITEN the data with the specified frequency
   * Fit data with FOURIER series, while frequency remains fixed. For frequency adjustment use refine, Fourier fit is prone to correlation of coefficients and frequency errors. 
   * RESET the plot scales to initial or last freeze by clicking button
   * Select periodogram method by ticking AMH|AOV|ATR|POW|WIN marks
   * Drag the SLIDER or type NH2 to modify no. of bins|harmonics*2
   * SAVE the (prewhitened?) DATA to a file using the File>Save Data menu
   * SAVE the PLOT to a file using the File>Save Plot menu
   * Use File>FREEZE & OPEN to freeze periodogram parameters and input new data 
   * Click on a bar to receive an informative message
        """
        QMessageBox.about(self, "About AOV Search", msg.strip())

    def on_author(self):  # Action on About>Help
        msg = """ AOV Period Search of Light Curve
        
        Author: Alex Schwarzenberg-Czerny, alex@camk.edu.pl
        Reference:
        Alex Schwarzenberg-Czerny, 1998, Baltic Astronomy, v.7, p.43-69.

        Licence: This python code is in the public domain. The associated Fortran routines are copyrighted (C) by Alex Schwarzenberg Czerny and put into your free use, provided they remain unchanged, and are not used for purposes carrying risk of any damages.
        History: 
        Origin: 19.10.2012
        Modiffied: 26.10.2012 - general streamline of windows
        Modiffied: 12.11.2012 - frequency error estimate in REFINE and NOTE box
        

        Acknowledgements:
        This routine embeds a matplotlib plot into PyQt4 GUI application following a demo by Eli Bendersky.
        Pyaov.py wrapper of f95 periodogram follows an example by Ewald Zietsman
        """
        QMessageBox.about(self, "AOV Search Package: History",
                          msg.strip())

        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
# =============Periodogram plot editing=============================
    def ax_loc(self, event, bound):
        low, up = getp(event.inaxes, bound)
        data = (event.xdata)
        if bound == 'ybound':
            data = (event.ydata)
        span = (up-low)
        frac = (data-low)/span
        if abs(frac) < self.marg:
            frac = 0.
        if abs(frac-1.) < self.marg:
            frac = 1.
        return low, span, frac

    def on_click(self, event):
        ax = event.inaxes
        if ax != None:
            yt = (self.ax_loc(event, 'ybound'))
            xt = (self.ax_loc(event, 'xbound'))
        else:
            return
        if event.button is MouseButton.LEFT:
            if self.last_click != None:
                print("Last RIGHT click ignored")
                self.last_click != None
            self.pick_point(event, yt, xt)
            return
        if event.button is MouseButton.RIGHT:
            ax.plot(xt[0]+xt[1]*xt[2], yt[0]+yt[1]*yt[2], 'r-o')
            if self.last_click != None:
                if ax is self.last_click[0].inaxes:
                    top = self.last_click[1][2] > 0.5
                    if ax is self.per.ax:
                        self.per.xr[0], self.per.xr[1] = self.pick_range(
                            event, self.last_click)
                    else:
                        self.lc.xr = self.pick_range(
                            event, self.last_click), 0.
                    self.on_draw(ax)
                    self.last_click = None
                else:
                    print(" no RIGHT clicks match, previous one ignored")
                    self.last_click = (event, yt, xt)
            else:
                self.last_click = (event, yt, xt)
        self.canvas.draw()

    def pick_point(self, event, yt, xt):
        if event.inaxes is self.lc.ax:
            # find nearby point
            xmarg = self.marg*xt[1]*0.1
            ymarg = self.marg*yt[1]
            b = np.nonzero(abs(self.ph.xd-event.xdata) < (2*xmarg))[0]
            r = np.sqrt(np.square((self.ph.xd[b] - event.xdata)/xmarg) +
                        np.square((self.ph.yd[b]-event.ydata)/ymarg))

            ind = r.argmin()
            if r[ind] < 1.:
                ind = self.phidx[b[ind]]
                self.lcd.ed[ind] = -self.lcd.ed[ind]
        else:
            n = len(self.per.xd)
            ib = int(n*0.035)+1
            ia = max(int((event.xdata-self.per.xr[0])/self.per.xr[2])-ib, 0)
            ib = min(ia+ib+ib, n-1)
            xm, fm, dx = pyaov.peak(self.per.yd[ia:ib])
            self.per.pxy = np.array(
                [xm*self.per.xr[2]+self.per.xd[ia], 0.5*fm])
            self.per.pplt = True
            self.per.ax.clear()
            self.per.ax.plot(self.per.xd, self.per.yd, linestyle='solid')
            if self.per.pplt:
                self.per.ax.plot(self.per.pxy[0], self.per.pxy[1], 'r-o')
            self.textbox.setText("%18.10g" % self.per.pxy[0])
        self.on_draw_lc()

    def pick_range(self, event, last_click):
        top = last_click[1][2] > 0.5
        xsav, ysav = (last_click[0].xdata), (last_click[0].ydata)
# zoom in periodogram between marks or...
        xup = max(xsav, event.xdata)
        xlow = min(xsav, event.xdata)
        if top:
            # ...zoom out periodogram by the span ratio
            xlow_old = (last_click[2][0])
            xspan = (last_click[2][1])
            xup_old = xlow+xspan
            ratio = xspan/abs(event.xdata-xsav)
            xup = xup+ratio*(xup_old-xup)
            xlow = max(xlow-ratio*(xlow-xlow_old), 0.)

        if event.inaxes is self.per.ax:
            self.per.xr[2] = min((xup-xlow)/300., self.per0.xr[2])
        return xlow, xup

    def on_reset(self):  # Reset intervalper & lc grids
        self.per.xr = np.array(self.per0.xr)
        self.lc.xr = np.array(self.lc0.xr)
        self.on_draw_per(self.radioGroup.checkedId())

    def on_refine(self):
        print('Refine\n')
        s = self.textbox.text()
        try:
            frout = float(s)
        except:
            print('Need number as note')
            return 0
        fr, dfr, self.lcd.yd, nh2 = \
            pyaov.refine(self.lcd.xd, self.lcd.yd, self.lcd.ed, frout)
        self.textbox.setText("%18.10g" % fr)
        print('refined with nh2= ', nh2)
        self.on_draw_per(self.radioGroup.checkedId())

    def on_prewhite(self):
        # prewhite data with fixed frequency by subtraction of
        # trig orthogonal polynomial series

        # frequency from box text
        print('Prewhite\n')
       # str = unicode(self.textbox.text())
        str = self.textbox.text()
        nu = float(str)
        if nu == 0.:
            print('on_prewhite: error-zero nu')

        fr, dfr, self.lcd.yd = pyaov.prew(self.lcd.xd, self.lcd.yd,
                                          self.lcd.ed, -nu, nh2=self.nh2)  # nu<0 not adjusted
        self.lcd.yd, self.lcd.nact, self.lcd.sw, self.lcd.av, self.lcd.vr = \
            pyaov.normalize(self.lcd.yd, self.lcd.ed, var=0.)

        # get periodogram
        self.on_draw_per(self.radioGroup.checkedId())

    def on_fourier(self):
        # prewhite data with fixed frequency by subtraction of
        # Fourier series

        # frequency from box text
        print('Fourier\n')

       # str = unicode(self.textbox.text())
        str = self.textbox.text()
        print(type(str), len(str), str)
        nu = float(str)
        if nu == 0.:
            print('on_fourier: error-zero nu\n')

        fr, dfr, self.lcd.yd, cof, dcof = pyaov.fouw(self.lcd.xd, self.lcd.yd,
                                                     self.lcd.ed, -nu, nh2=self.nh2)  # nu<0 not adjusted
        self.lcd.yd, self.lcd.nact, self.lcd.sw, self.lcd.av, self.lcd.vr = \
            pyaov.normalize(self.lcd.yd, self.lcd.ed, var=0.)

        # get periodogram
        self.on_draw_per(self.radioGroup.checkedId())

    def on_draw(self, ax):
        if ax is self.per.ax:
            self.on_draw_per(self.radioGroup.checkedId())
        if ax is self.lc.ax:
            self.on_draw_lc()

    def on_draw_lc(self):
        """ Redraws the folded light curve
        """
        # data from box text
       # str = unicode(self.textbox.text())
        str = self.textbox.text()
        lcd = self.lcd
        nu = float(str)
        if nu == 0.:  # print raw light curve
            nu = 1./(lcd.xd.max()-lcd.xd.min())

        # clear the axes and redraw the plot anew
        self.lc.ax.clear()
        self.lc.ax.grid(self.grid_cb.isChecked())
        imin = lcd.yd.argmin()
        x = (lcd.xd-lcd.xd[imin])*nu+0.5
        x = x-np.floor(x)
        y = np.array(lcd.yd)
        e = np.array(lcd.ed)
        if (self.lc.xr[1] > 1.):
            x = np.append(x, np.add(x, 1.))
            y = np.append(y, y)
            e = np.append(e, e)
        # dl. 1.20, wart. 0<2len
        idx = np.nonzero(np.logical_and(
            x > self.lc.xr[0], x < self.lc.xr[1]))[0]

        ph = My_Lc()
        ph.xd = np.array(x[idx])
        ph.yd = np.array(y[idx])
        ph.ed = np.array(e[idx])
        idx = idx[idx % lcd.xd.size] % lcd.xd.size

        b = ph.ed > 0.
        self.lc.ax.plot(ph.xd[b], ph.yd[b], 'g o')
        b = ph.ed <= 0.
        self.lc.ax.plot(ph.xd[b], ph.yd[b], 'r o')
        self.lc.ax.set_xlim(self.lc.xr[0], self.lc.xr[1])
        yup = ph.yd.max()
        ylow = ph.yd.min()
        yspan = (yup-ylow)*(0.5+3*self.marg)
        yav = (yup+ylow)*0.5
        self.lc.ax.set_ylim(yav-yspan, yav+yspan)

        self.canvas.draw()
        del x, y, e
        self.ph, self.phidx = ph, idx

    def on_slide_nh2(self):
        nh2 = max(int(10.**(self.slider.value()/100.*1.6)+0.5)+3, 3)
        self.nh2box.setText("%3d" % nh2)
        self.canvas.draw()

    def on_edit_nh2(self):
        #str = unicode(self.nh2box.text())
        s = self.nh2box.text()
        nh2 = max(int(s), 3)
        if (nh2 != self.nh2):
            self.nh2 = (nh2)
            self.nh2box.setText("%3d" % nh2)
            self.on_draw_per(self.radioGroup.checkedId())

    def on_draw_per(self, meth):
        """ Redraws the periodogramme
        """
#        meth = self.radioGroup.checkedId()
        #str = unicode(self.nh2box.text())
        s = self.nh2box.text()
        self.nh2 = max(int(s), 3)
        self.nh2box.setText("%3d" % self.nh2)
        print('method= ', 'amh aov atr pow win'.split()[meth], '\n')

        self.per.ax.clear()
        self.per.ax.grid(self.grid_cb.isChecked())
        # Calculate perioddogram
        if meth == 0:
            self.per.yd, self.per.xd, self.per.pxy[0] = \
                pyaov.amhw(self.lcd.xd, self.lcd.yd, self.lcd.ed,
                           self.per.xr[1], self.per.xr[2], fr0=self.per.xr[0], nh2=self.nh2)
        elif meth == 1:
            self.per.yd, self.per.xd, self.per.pxy[0] = \
                pyaov.aovw(self.lcd.xd, self.lcd.yd, self.lcd.ed,
                           self.per.xr[1], self.per.xr[2], fr0=self.per.xr[0], nh2=self.nh2)
        elif meth == 2:
            self.per.yd, self.per.xd, self.per.pxy[0] = \
                pyaov.atrw(self.lcd.xd, self.lcd.yd, self.lcd.ed,
                           self.per.xr[1], self.per.xr[2], fr0=self.per.xr[0], nh2=self.nh2)
        elif meth == 3:
            self.per.yd, self.per.xd, self.per.pxy[0] = \
                pyaov.pspw(self.lcd.xd, self.lcd.yd, self.lcd.ed,
                           self.per.xr[1], self.per.xr[2], fr0=self.per.xr[0])
        elif meth == 4:
            self.per.yd, self.per.xd, self.per.pxy[0] = \
                pyaov.pspw(self.lcd.xd, np.ones(len(self.xd), dtype=float),
                           self.lcd.ed, self.per.xr[1], self.per.xr[2], fr0=self.per.xr[0])
        else:
            print('on_draw_per: Wrong mode')

        # set new frequency and recalculate lc
        self.per.pxy[1] = 0.5*np.amax(self.per.yd)
        self.per.pplt = True
        self.per.ax.plot(self.per.xd, self.per.yd, linestyle='solid')
        if self.per.pplt:
            self.per.ax.plot(self.per.pxy[0], self.per.pxy[1], 'r-o')
        self.canvas.draw()
        self.textbox.setText("%18.10g" % self.per.pxy[0])
        self.on_draw_lc()

    def on_note(self):
        #str = unicode(self.notebox.text())
        s = self.notebox.text()
        if (s != self.strold and s != " "):
            print(s)
            self.notebox.setText(" ")
            self.canvas.draw()
            self.strold = s

    def create_main_frame(self):  # Organize the main window
        self.main_frame = QWidget()

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        screen = QDesktopWidget().screenGeometry()
        #size = self.geometry()
        self.fig = Figure((screen.width()*0.8/self.dpi,
                           screen.height()*0.5/self.dpi), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.per0.ax = self.fig.add_subplot(211)
        self.lc0.ax = self.fig.add_subplot(212)

        # Bind the 'pick' event for clicking on one of the bars

        cid = self.fig.canvas.mpl_connect('button_press_event',
                                          self.on_click)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

# Row of GUI view controls ===========================

        # edited frequency field
        freq_label = QLabel("freq=")
        self.textbox = QLineEdit()
        self.textbox.setMinimumWidth(200)
        self.textbox.editingFinished.connect(self.on_draw_lc)
      #  self.connect(self.textbox, SIGNAL('editingFinished ()'), \
       #          self.on_draw_lc)

        # refine frequency button
        self.refine_button = QPushButton("&Refine")
        self.refine_button.clicked.connect(self.on_refine)
        # self.connect(self.refine_button, SIGNAL('clicked()'), \
        #   self.on_refine)

        # prewhite button
        self.prewhite_button = QPushButton("&Prewhite")
        self.prewhite_button.clicked.connect(self.on_prewhite)
       # self.connect(self.prewhite_button, SIGNAL('clicked()'), \
        #    self.on_prewhite)

        # Fourier coefficients button
        self.fourier_button = QPushButton("&Fourier")
        self.fourier_button.clicked.connect(self.on_fourier)
      #  self.connect(self.fourier_button, SIGNAL('clicked()'), \
       #     self.on_fourier)

        # zoom out button
        self.zoom_button = QPushButton("&Reset")
        self.zoom_button.clicked.connect(self.on_reset)
     #   self.connect(self.zoom_button, SIGNAL('clicked()'), self.on_reset)

        # show grid button
        self.grid_cb = QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.grid_cb.stateChanged.connect(self.on_draw_per)
      #  self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), \
       #     self.on_draw_per)

        # Layout of GUI view controls bar
        hbox = QHBoxLayout()

        for w in [freq_label, self.textbox, self.refine_button,
                  self.prewhite_button, self.fourier_button,
                  self.zoom_button, self.grid_cb]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox)

# Row of method controls as RadioButtons =========================
        listOfChoices = 'amh aov atr pow win'.split()
        self.radioGroup = QButtonGroup()
        self.radioGroup.setExclusive(True)
        rbox = QHBoxLayout()

        bIsFirst = True
        nMaxLen = 0
        for i, row in enumerate(listOfChoices):
            radio = QRadioButton(row)
            self.radioGroup.addButton(radio, i)
            if bIsFirst:
                radio.setChecked(True)
                bIsFirst = False
            if len(row) > nMaxLen:
                nMaxLen = len(row)
# self.listWidget.itemDoubleClicked.connect(self.showItem) <<< new syntax

            radio.toggled.connect(self.on_draw_per)
         #  self.connect(radio, SIGNAL("toggled(bool)"), \
            #       self.on_draw_per)
            rbox.addWidget(radio)
            rbox.setAlignment(radio, Qt.AlignVCenter)

        # slider field
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 100)
        self.slider.setValue(0)
        self.slider.setTracking(True)  # may call while mooving
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.slider.valueChanged.connect(self.on_slide_nh2)
       # self.connect(self.slider, SIGNAL('valueChanged(int)'), \
        #     self.on_slide_nh2)
        rbox.addWidget(self.slider)
        rbox.setAlignment(self.slider, Qt.AlignVCenter)

        # edited nh2 field
        nh2_label = QLabel("nh2=")
        rbox.addWidget(nh2_label)
        rbox.setAlignment(nh2_label, Qt.AlignVCenter)
        self.nh2box = QLineEdit()
        self.nh2box.setMinimumWidth(30)
        self.nh2box.editingFinished.connect(self.on_edit_nh2)
       # self.connect(self.nh2box, SIGNAL('editingFinished ()'), \
        #         self.on_edit_nh2)
        rbox.addWidget(self.nh2box)
        rbox.setAlignment(self.nh2box, Qt.AlignVCenter)

        vbox.addLayout(rbox)

# Row for notes  ============================================
        note_label = QLabel("Note:")
        self.notebox = QLineEdit()
        # self.notebox.setMinimumWidth(screen.width()*0.8)
        self.notebox.setMinimumWidth(int(screen.width()*0.8))  # KPI
        self.notebox.editingFinished.connect(self.on_note)
       # self.connect(self.notebox, SIGNAL('editingFinished ()'), \
        #         self.on_note)

        nbox = QHBoxLayout()
        nbox.addWidget(note_label)
        nbox.setAlignment(note_label, Qt.AlignVCenter)
        nbox.addWidget(self.notebox)
        nbox.setAlignment(self.notebox, Qt.AlignVCenter)

        vbox.addLayout(nbox)

# End menu rows ==========================================

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("by Alex Schwarzenberg-Czerny, 2012")
        self.statusBar().addWidget(self.status_text, 1)

    def create_menu(self):  # Create File> and Help> Menus
        self.file_menu = self.menuBar().addMenu("&File")
        open_file_action = self.create_action("&Open file",
                                              shortcut="Ctrl+O", slot=self.open_dialog,
                                              tip="Open ascii data file (time, value columns)")

        freeze_file_action = self.create_action("&Freeze and Open",
                                                shortcut="Ctrl+F", slot=self.freeze_dialog,
                                                tip="Freeze parameters & open ascii data file")
        out_data_action = self.create_action("&Save data",
                                             shortcut="Ctrl+D", slot=self.save_data,
                                             tip="Save the prewhitened(?) data")
        out_period_action = self.create_action("&Save periodogram",
                                               shortcut="Ctrl+P", slot=self.save_period,
                                               tip="Save the current periodogram values")
        out_fig_action = self.create_action("&Save plot",
                                            shortcut="Ctrl+S", slot=self.save_plot,
                                            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close,
                                         shortcut="Ctrl+Q", tip="Close the application")

        self.add_actions(self.file_menu,
                         (open_file_action, freeze_file_action, out_data_action,
                          out_period_action, out_fig_action, None, quit_action))

        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",
                                          shortcut='F1', slot=self.on_about,
                                          tip='About the demo')
        author_action = self.create_action("A&uthor",
                                           shortcut='Ctrl+A', slot=self.on_author,
                                           tip='Author & History')
        self.add_actions(self.help_menu, (about_action, author_action))

    def add_actions(self, target, actions):  # Menus support
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None,
                      icon=None, tip=None, checkable=False,
                      signal="triggered()"):
        # more Menus support
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            #self.connect(action, SIGNAL(signal), slot)
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action

    def input_lc(self, filename):  # input data from file
        # read data file
        #z = atpy.Table.read(filename, type='ascii', Reader=asc.NoHeader)
        extension = os.path.splitext(filename)[1][1:]

        if extension in ['csv', 'dat', 'txt', 'asc', 'ascii']:
            data = ascii.read(filename)
        elif extension in ['fits']:
            data = fits.getdata(filename)
        else:
            raise IOError("Invalid extension for filename=", filename)

        valid_x = ['col1', 'mjd', 'MJD', 'hjd', 'HJD']
        valid_y = ['col2', 'Flux', 'flux', 'FLUX', 'MAG', 'mag', 'Mag']
        valid_e = ['col3', 'Flux_err', 'Fluxerr', 'flux_err', 'fluxerr',
                   'FLUX_ERR', 'FLUXERR', 'MAG_ERR', 'MAGERR', 'mag_err', 'magerr',
                   'Mag_err', 'Magerr']

        try:
            lc = My_Lc()
            lc.fname = filename
            xcol = next(i for i in valid_x if i in data.columns)
            ycol = next(i for i in valid_y if i in data.columns)
            lc.xd = np.array(data[xcol])
            lc.yd = np.array(data[ycol])
        except:
            raise ValueError("No valid columns names in ", data.columns)
        try:
            ecol = next(i for i in valid_e if i in data.columns)
            lc.ed = np.array(data[ecol])
        except:
            lc.ed = np.ones(lc.xd.size, dtype=float)
        del data
        return lc


def main():  # Main execution loop
    #    sys.exit()
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
        return lc


def main():  # Main execution loop
    #    sys.exit()
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
