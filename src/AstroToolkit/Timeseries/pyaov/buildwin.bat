:: This is a windows batch file for building pyaov and _aov
:: Intel fortran is a prerequisite
:: modified keith inight 6/7/24 to use wrapping rather than f2py

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

ifx -Qopenmp aovconst.f90 aovsub.f90 aov.f90 aovwrapper.f90 -o aovwrapper -fpp
path > path.dat
IF %ERRORLEVEL% NEQ 0 (
	Echo:
	Echo =======================================
	Echo Compilation Failed. An error was found.
	Echo =======================================
	Echo:
	EXIT
	)
	Echo:
	Echo =======================
	Echo Compilation Successful.
	Echo =======================
	Echo:

del *.mod
del *.obj
