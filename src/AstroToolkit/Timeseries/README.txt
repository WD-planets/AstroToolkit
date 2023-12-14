AOVGUI Python Package
**** conversion to Python3 unfinished

Purpose: For time and frequency space analysis of unevenly-sampled time series.
==============================================================================
The package consists of two parts: pyaov.py supplies individual python routines intended for advanced work, aovgui.py provides their gui frontend for interactive work and is intended for beginners.

(C) Alex Schwarzenberg-Czerny 2011-2012 alex@camk.edu.pl
The whole package aov is copyrighted. It may be distributed free,
provided that no changes are made to any file and its copyright statements.
This is so as to prevent unmanagable forking of the package.

Requirements
====================================================
Both packages require (free) gfortran & f2py fortran and fortran-to-python compiler and python 2.7 installation, i.e. the latest stable one. They also need python matplotlib.pyplot and numpy packages, usually installed by default. Additionally, aovgui.py requires other popular python packages: sys, os, random, scipy and pylab and also less popular ones: atpy, asciitable, pprint. I am not expert on python installation so no detailed instruction here, though you may have tu use pip or similar python installer. Python widget library PyQt4 has to be installed via linux, perhaps by 
>sudo apt-get install python-qt4 qt4-dev-tools python-qt4-dev pyqt4-dev-tools
or similar action, same appliues to f2py compiler/wrapper, if it is not already installed.

Note: at present the AOVGUI package is configured for single precision
observations(magnitudes, velocities etc.). To change it to double precision
pls edit aovconst.f90 so that SP=8 and CP=16.

Installation
==========================================================
After untarring the distribution file you may need to enable execution
of the builder script: chmod a+x build.sh. Next execute it by issuing from shell: ./build.sh, with long output.

Now you are ready for a test run of pyaov.py and aovgui.py,
directly from python. You need to imort them or else signal to python 
they are needed.

aovgui.py package
====================================================
Purpose:
GUI for interactive analysis of unevenly sampled time series

Test run: Start aovgui.py by typing from the shell
python aovgui.py | tee aov.log
This way you will get copy of screen messages into a aov.log file.
Select your data file (or ex1.dat) from the popping out window
and read a brief instruction from the Help>About. 

Additionally, source of aovgui.py may serve an (admittedly complex)
example how to call the pyaov.py routines in python scripts.

For detailed requirements see the import/from statements in the aovgui.py header.

pyaov.py routines
==================================================== 
Purpose:
Routines for analysis of unevenly sampled time series

Test run: To start an example type from the shell: python ex1.py
To proceed you have to close each subsequent graphic window (press cross in its 
right upper corner). The same applies for ex2.py

Now you are set to use PyAOV package. Enter interactive python mode and 
examples to proceed in a similar way. Use online help as indicated 
in pyaov.py header and consult documentation from the latex pyaov.tex file


