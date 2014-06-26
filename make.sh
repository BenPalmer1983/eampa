#!/bin/bash
mpif90 -o eampa.x src/kinds.f90 src/constants.f90 src/mpif.f90 src/strings.f90 src/maths.f90 src/general.f90 src/units.f90 src/data.f90 src/initialise.f90 src/input.f90 src/prep.f90 src/output.f90 src/calc.f90 src/optimise.f90 src/prepeam.f90 src/prepdl.f90 src/pwbatch.f90 src/run.f90 src/eampa.f90 
sleep 1 
rm *.mod