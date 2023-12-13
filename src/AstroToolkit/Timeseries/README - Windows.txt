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
Both packages require (free) Intel fortran, Visual Studio and a python installation. They also need python matplotlib.pyplot and numpy packages, usually installed by default. Additionally, aovgui.py requires other popular python packages: sys, os, random, scipy and pylab and also less popular ones: atpy, asciitable, pprint,QT5. All of these are automatically included in the Anaconda distribution.  

Note: at present the AOVGUI package is configured for single precision
observations(magnitudes, velocities etc.). To change it to double precision
pls edit aovconst.f90 so that SP=8 and CP=16.

Installation
==========================================================
Create a folder for the code.
Next untar the distribution file and copy the contents to the folder. 
Then open a command prompt and change the default directory to this folder using CD {folder}.
Then type "build.bat"
	A long console list of output from f2py and the fortran and C compilers is then displayed. This includes many apparent errors which are most likley not a problem. 
	If the run is successful a message "f2py was successful" is displayed and a test program is run that displays a number of periodograms.


The intel fortran environment requires a number of environment variables to be set in order to work. Intel have created setvars.bat to automate this. Setvars.bat needs to be run before using aovgui or any other python program accessing the aov package. This can be automated by firstly creating a shortcut to setvars.bat (which resides in C:\Program Files (x86)\Intel\oneAPI if the default installation was used). This should then be copied to the startup folder which can be found by pressing the windows button and "R"  and then typing placing it in either    "shell:common startup" without quotes. 
	

Now you are ready for a test run of pyaov.py and aovgui.py,
directly from python. You need to imort them or else signal to python 
they are needed.

aovgui.py package
====================================================
Purpose:
GUI for interactive analysis of unevenly sampled time series

Test run: Start aovgui.py by typing from the shell
python aovgui.py 

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


