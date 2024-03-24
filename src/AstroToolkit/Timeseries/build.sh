#!/bin/bash

# script to build pyaov and _aov
# alex schwarzenberg-czerny alex@camk.edu.pl 2012

# to run: sh build.sh

# usage: run command1 && run command2 && run command3

run() {
  $@
  if [ $? -ne 0 ]
  then
    tput setaf 2; echo "$* failed with exit code $?"; tput sgr0
    echo 
    return 1
  else
    tput setaf 2; echo "finished OK"; tput sgr0
    return 0
  fi
}

# first clean old executables & libraries
rm *.pyc *.so *.mod *.toc *.aux *.log

# next compile
run python -m numpy.f2py -c --f90flags=-fopenmp -lgomp aovconst.f90 aovsub.f90 aov.f90 -m aov
# finished
