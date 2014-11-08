#!/bin/bash
mkdir -p bin
module load /gpfs/apps/modules/apps/openmpi/v1.6.5/gcc-tm-ib/v4.8.2
mpif90 -g -fcheck=all -Wall -mtune=native -o bin/eampa.x \
src/kinds.f90 src/msubs.f90 src/constants.f90 \
src/maths.f90 src/general.f90 src/units.f90 src/globals.f90  \
src/initialise.f90 src/loadData.f90 src/output.f90 \
src/readInput.f90 src/readEAM.f90 src/makeConfig.f90 src/readConfig.f90  \
src/neighbourList.f90 src/prepEAM.f90 src/prep.f90  \
src/calcEAM.f90  \
src/relax.f90  \
src/bulkProperties.f90 src/testEAM.f90 src/calcEval.f90  \
src/optimise.f90  \
src/pwBatch.f90 src/clean.f90  \
src/eampa.f90
sleep 1 
rm *.mod