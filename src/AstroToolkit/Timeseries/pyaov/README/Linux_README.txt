Purpose
=======
For time and frequency space analysis of unevenly-sampled time series.



Note
====
The package consists of  pyaov.py which supplies individual python methods by accessing compiled fortran subroutines.

(C) Alex Schwarzenberg-Czerny 2011-2012 alex@camk.edu.pl
The whole package aov is copyrighted. It may be distributed free,
provided that no changes are made to any file and its copyright statements.
This is so as to prevent unmanagable forking of the package.

Modified by Keith Inight 6/7/24 to use wrapping rather than f2py as the fortran interface.
Modified by Ethan Moorfield for integration with the AstroToolkit Python package



Requirements
============
The package requires (free) GCC gfortran.



Installation
============
This guide assumes you already have Python installed, with the matplotlib and numpy packages. These are included in ATK.

1. Install GCC gfortran. In most cases, this should be downloaded from https://gfortran.meteodat.ch/download/

2. Run the following command from a terminal: python -c "from AstroToolkit.Setup import tsbuild; tsbuild()"

    If compilation is successful, pyaov's installation is complete.