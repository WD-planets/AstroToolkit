#!/bin/bash

# script to build aov routines
# alex schwarzenberg-czerny alex@camk.edu.pl 2012
# modified by keith inight 6/7/24 to use wrapping rather than f2py
# modified by Ethan Moorfield for integration with the AstroToolkit package

run() {
  $@
  f95 -fopenmp -lgomp aovconst.f90 aovsub.f90 aov.f90 aovwrapper.f90 -o aovwrapper -static -cpp
  if [ $? != 0 ]
  then
    echo 
	echo =======================================
	echo Compilation Failed. An error was found.
	echo =======================================
	echo 
  else
    echo 
    echo =======================
	echo Compilation Successful.
	echo =======================
	echo 
  fi
}

run

