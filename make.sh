#!/bin/bash
echo "Current working directory:"
pwd
echo "Make bin directory"
mkdir -p bin
echo "Check formatting/layout of source files"
#python make.py 1
echo "Compile source"
mkdir -p bin
#mpif90 -g -fcheck=all -Wall -mtune=native -std=f95 -o
mpif90 -g -fcheck=all -Wall -mtune=native -Wunused-function -o bin/eampa.x \
src/kinds.f90 src/msubs.f90 src/constants.f90 \
src/maths.f90 src/mMaths.f90 src/general.f90 src/units.f90 src/globals.f90  \
src/initialise.f90 src/loadData.f90 src/output.f90 \
src/readInput.f90 src/readEAM.f90 src/readConfig.f90 \
src/neighbourList.f90 src/preCalc.f90 src/calcEAM.f90 \
src/finalise.f90 \
src/eampa.f90
#src/eampa.f90 &> logs/make.log
#sleep 1 
rm *.mod
#echo "Check compile result"
#python make.py 2