#!/bin/bash
#Make the binary output directory
echo "Make bin directory"
mkdir -p bin
#Check source files are OK, remove tabs
echo "Check formatting/layout of source files"
python make.py 1
sleep 1 
#Compile first time  mpif90 -O3 -o 
echo "Compile code - attempt 1"
#mpif90 -g -fcheck=all -Wall -mtune=native -std=f95 -o
mpif90 -g -fcheck=all -Wall -mtune=native -o bin/eampa.x \
src/kinds.f90 src/msubs.f90 src/constants.f90 \
src/maths.f90 src/mMaths.f90 src/general.f90 src/units.f90 src/globals.f90  \
src/initialise.f90 src/loadData.f90 src/output.f90 \
src/readInput.f90 src/readEAM.f90 src/makeConfig.f90 src/readConfig.f90  \
src/neighbourList.f90 src/prepEAM.f90 src/prep.f90  \
src/calcEAM.f90 \
src/relax.f90 \
src/bulkProperties.f90 src/testEAM.f90  \
src/calcEval.f90  \
src/optimise.f90  \
src/pwBatch.f90 src/clean.f90  \
src/eampa.f90 &> logs/build.log
sleep 1 
rm *.mod
#Remove unused variables
echo "Remove unused variables from source files"
python make.py 3
echo "Check formatting/layout of source files"
python make.py 1
#Recompile
mpif90 -g -fcheck=all -Wall -mtune=native -o bin/eampa.x \
src/kinds.f90 src/msubs.f90 src/constants.f90 \
src/maths.f90 src/mMaths.f90 src/general.f90 src/units.f90 src/globals.f90  \
src/initialise.f90 src/loadData.f90 src/output.f90 \
src/readInput.f90 src/readEAM.f90 src/makeConfig.f90 src/readConfig.f90  \
src/neighbourList.f90 src/prepEAM.f90 src/prep.f90  \
src/calcEAM.f90 \
src/relax.f90 \
src/bulkProperties.f90 src/testEAM.f90  \
src/calcEval.f90  \
src/optimise.f90  \
src/pwBatch.f90 src/clean.f90  \
src/eampa.f90 &> logs/make.log
sleep 1 
rm *.mod 
echo "Check compile result"
python make.py 2