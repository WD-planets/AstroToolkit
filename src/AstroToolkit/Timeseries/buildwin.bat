
rem This is a windows batch file for building pyaov and _aov
rem Intel fortran and visual studio C are prerequisites

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

rem -c builds an extension module

python -m numpy.f2py -c --f90flags=-Qopenmp aovconst.f90 aovsub.f90 aov.f90 -m aov
IF %ERRORLEVEL% NEQ 0 (
	Echo
	Echo ********************************
	Echo
	Echo f2py failed An error was found
	Echo
	Echo ********************************
	)
	Echo ********************************
	Echo
	Echo f2py was successful
	Echo
	Echo ********************************
setx KMP_DUPLICATE_LIB_OK TRUE
