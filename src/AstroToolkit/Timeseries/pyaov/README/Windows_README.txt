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
The package requires (free) Intel fortran and Microsoft Visual Studio Build Tools. 



Installation
============
This guide assumes you already have python installed, with the matplotlib and numpy packages. These are included in ATK. All following install locations should be left to their defaults.

1. Install Microsoft Build Tools
   These can be most easily installed by going to https://visualstudio.microsoft.com/downloads/ and searching for "Build Tools for Visual Studio 2022"
   - Before installing, you must make sure to include "Desktop Development with C++" under workloads. All optional modules can be disabled except for MSVC and the Windows SDK.

2. Install Intel Fortran
   This can be installed from https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran
   A warning will state that this will only work from the command line, this is fine.

3. Run the following command in a terminal: python -c "from AstroToolkit.Setup import tsbuild; tsbuild()".
	
   If compilation is successful, pyaov's installation is complete.



Known Issues
============

1. Microsoft Build Tools not found
   If compilation fails due to not being able to find your build tools, you must add "VS2022INSTALLDIR" to your system's environment variables with the
   path to them: "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools".